# This file is part of Bertini 2.
#
# python/bertini/config.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/config.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/bertini/config.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""
Reusable, low-friction ergonomics for Bertini config structs.

The C++ bindings expose each config struct (``SteppingConfig``, ``NewtonConfig``,
the endgame configs, the zero-dim configs, ...) with one writable attribute per
field, and expose each *owner* (every tracker and nag_algorithm) with a uniform,
type-list-driven config interface: ``set_config``, ``get_config`` and
``config_types`` (see ``configured_visitor.hpp``).

This module layers Pythonic conveniences on top of those, with **no per-config
and no per-algorithm code** -- everything is discovered by introspection, so a
newly added algorithm and its configs get the whole interface for free:

  * configs gain ``update(**kwargs)`` (chainable, validated), ``to_dict()``,
    ``from_dict()``, a readable ``repr`` and value-equality;
  * owners gain ``configure(**kwargs)`` and ``config_names()``, built on the
    owner's own ``config_types()`` so they always reflect what the owner accepts.
"""

import copyreg
import re


# ---------------------------------------------------------------------------
# field discovery
# ---------------------------------------------------------------------------

def writable_fields(cls):
    """The names of the writable (def_readwrite) data members of a bound class.

    Discovered by walking the class for data descriptors (``__set__``), so it
    requires no hand-maintained field lists and tracks the C++ struct exactly.
    """
    fields = []
    for klass in getattr(cls, "__mro__", (cls,)):
        for name, val in vars(klass).items():
            if name.startswith("_") or name in fields:
                continue
            descr = type(val)
            if not hasattr(descr, "__set__") or not hasattr(descr, "__get__"):
                continue
            # a read-only python property would have __set__ on its type but no fset
            if isinstance(val, property) and val.fset is None:
                continue
            fields.append(name)
    return tuple(fields)


def _camel_to_snake(name):
    s = re.sub(r"(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", s).lower()


def config_key(cls):
    """A short, snake_case keyword for a config class (drops a trailing 'Config').

    e.g. SteppingConfig -> 'stepping', AMPConfig -> 'amp',
         PostProcessingConfig -> 'post_processing'.
    """
    name = cls.__name__
    if name.endswith("Config"):
        name = name[: -len("Config")]
    return _camel_to_snake(name)


# ---------------------------------------------------------------------------
# config-class enhancement
# ---------------------------------------------------------------------------

def _coerced_setattr(obj, key, value):
    """setattr, retrying a string as the field's numeric type.

    Strings are the blessed noise-free way to express numeric settings
    (e.g. max_step_size="0.05", final_tolerance="1e-11"); plain Python floats
    stay rejected for the *exact* fields by policy, since 0.05 the double is not
    1/20.  The catch is that the fields have different underlying types:
    mpq_rational / mpfr_float fields (step sizes, factors) take a multiprecision
    Float, while NumErrorT / double fields (tolerances, AMP bounds, min_step_size)
    take a plain double, and integer-count fields take an int.  So one string
    spelling works for EVERY field: we try it as an exact Float first, then as a
    plain float, then as an int, and keep the first that the field accepts.

    One family of fields is stored as an exact RATIONAL rather than as a float -- a sample
    ladder's ratio, a homotopy's times -- and those read only the fraction spelling, so
    ``_exact_to_rational`` is among the candidates: it takes ``'1/2'`` and ``'0.5'`` alike and
    makes the rational each one denotes.  An int is coerced for the same reason (a rational
    field refuses a bare Python int).

    Coercion only happens for exact spellings (a float that the field rejects raises straight
    through, with a message saying why), so the float-rejection policy on the exact fields
    stands: passing the double 0.05 to max_step_size is still an error.
    """
    from fractions import Fraction

    try:
        setattr(obj, key, value)
        return
    except TypeError as first_error:        # Boost.Python.ArgumentError is a TypeError
        if isinstance(value, float):
            raise TypeError(
                "{0!r} does not take the Python float {1!r}: it stores an exact value, and "
                "0.1 the double is not one tenth.  Pass an exact spelling instead -- a string "
                "('1/10' or '0.1'), an int, or a fractions.Fraction.".format(key, value)
            ) from first_error
        if not isinstance(value, (str, int, Fraction)):
            raise

    from .multiprec import real_mp
    from ._coefficients import _exact_to_rational

    # A string could be any of the field types, so try them in turn; an int or a Fraction only
    # ever reaches here because the field wanted an exact rational.
    candidates = (real_mp, _exact_to_rational, float, int) if isinstance(value, str) \
        else (_exact_to_rational,)

    last_error = None
    for convert in candidates:
        try:
            # int("1e-7") is a ValueError, real_mp("1/2") a RuntimeError out of the mpfr string
            # parser -- neither says anything about the field, so keep going rather than let an
            # incidental parse failure escape as the diagnosis
            converted = convert(value)
        except (ValueError, TypeError, RuntimeError):
            continue
        try:
            setattr(obj, key, converted)
            return
        except TypeError as e:              # field rejected this representation; try the next
            last_error = e
    raise last_error if last_error is not None else TypeError(
        "could not set {0!r} to {1!r}".format(key, value))


def _did_you_mean(name, candidates):
    """A ``  Did you mean 'x'?`` clause for a misspelling, or empty when nothing is close.

    A settings call names fields as keywords, so a typo is silent at the call site and only the
    error message can catch it.  Listing every valid name says what is available but not what
    the user probably meant, and these lists run to dozens of names.  Cutoff 0.6, stdlib
    difflib: close enough to catch a transposition or a missing underscore, far enough not to
    guess wildly.
    """
    from difflib import get_close_matches

    close = get_close_matches(str(name), [str(c) for c in candidates], n=3, cutoff=0.6)
    if not close:
        return ""
    if len(close) == 1:
        return "  Did you mean {0!r}?".format(close[0])
    return "  Did you mean one of {0}?".format(", ".join(repr(c) for c in close))


def _make_update(fields):
    def update(self, **kwargs):
        """Set one or more fields at once; returns self so calls can chain.

        Numeric fields accept strings (e.g. max_step_size="0.05"), which are
        converted exactly to multiprecision values.  Raises AttributeError on
        an unknown/misspelled field name.
        """
        for key, value in kwargs.items():
            if key not in fields:
                raise AttributeError(
                    "{0} has no config field {1!r}.{2}  Valid fields: {3}".format(
                        type(self).__name__, key, _did_you_mean(key, fields), list(fields)))
            _coerced_setattr(self, key, value)
        return self
    return update


def _make_to_dict(fields):
    def to_dict(self):
        """The config's fields as an ordinary dict of {name: value}."""
        return {k: getattr(self, k) for k in fields}
    return to_dict


def _make_from_dict():
    def from_dict(cls, mapping):
        """Build a config from a dict (default-constructs, then update())."""
        return cls().update(**dict(mapping))
    return classmethod(from_dict)


def _make_repr(fields):
    def __repr__(self):
        parts = []
        for k in fields:
            try:
                parts.append("{0}={1!r}".format(k, getattr(self, k)))
            except Exception:
                pass
        return "{0}({1})".format(type(self).__name__, ", ".join(parts))
    return __repr__


def _make_eq(fields):
    def __eq__(self, other):
        if type(self) is not type(other):
            return NotImplemented
        return all(getattr(self, k) == getattr(other, k) for k in fields)
    return __eq__


def _reconstruct_enhanced(cls, mapping):
    """Pickle reconstructor for enhanced config/metadata classes: rebuild from a field dict."""
    return cls.from_dict(mapping)


def _make_reduce():
    def __reduce__(self):
        # Pickle via the (to_dict / from_dict) round-trip.  Works for any enhanced class; the field
        # values must themselves be picklable (numbers, including multiprec.Complex/Float).
        return (_reconstruct_enhanced, (type(self), self.to_dict()))
    return __reduce__


def _make_copy():
    # copy.copy / copy.deepcopy rebuild through the same to_dict/from_dict round-trip as pickle.
    # Defining __deepcopy__ explicitly is load-bearing: the generic copy.deepcopy fallback would
    # deep-copy the to_dict() mapping's *values* (the multiprecision scalars), which aliases/corrupts
    # them on some interpreters (observed: Python 3.10 collapses the mpq_rational fields).  The field
    # values are immutable scalars re-fetched from the getters, so rebuilding from them -- without
    # deep-copying the scalars at all -- is both correct and portable.
    def __copy__(self):
        return type(self).from_dict(self.to_dict())

    def __deepcopy__(self, memo):
        new = type(self).from_dict(self.to_dict())
        memo[id(self)] = new
        return new

    return __copy__, __deepcopy__


def _enhance_config_class(cls):
    """Add update/to_dict/from_dict/repr/eq/pickle to a bound config class (idempotent)."""
    if getattr(cls, "_b2_config_enhanced", False):
        return cls
    fields = writable_fields(cls)
    cls._b2_fields = fields
    cls.update = _make_update(fields)
    cls.set = cls.update           # `set` reads more naturally to many users; same behavior
    cls.to_dict = _make_to_dict(fields)
    cls.from_dict = _make_from_dict()
    try:
        cls.__repr__ = _make_repr(fields)
        cls.__eq__ = _make_eq(fields)
        # Make the class picklable (copy.copy/deepcopy too) via to_dict/from_dict, overriding the
        # Boost.Python "pickling not enabled" default.
        cls.__reduce__ = _make_reduce()
        cls.__copy__, cls.__deepcopy__ = _make_copy()
    except (TypeError, AttributeError):
        # some bound types may refuse dunder assignment; the rest still apply
        pass
    cls._b2_config_enhanced = True
    return cls


# a class is an "owner" (tracker/algorithm), not a config, if it carries the
# type-list-driven config interface from ConfiguredVisitor.
def _is_owner(cls):
    return hasattr(cls, "set_config") and hasattr(cls, "config_types")


def _looks_like_config(cls):
    return (isinstance(cls, type)
            and not _is_owner(cls)
            and not _is_bound_enum(cls)
            and len(writable_fields(cls)) > 0)


def _is_bound_enum(cls):
    """True for a Boost.Python enum class.

    Boost.Python enums subclass ``int``, so ``writable_fields`` would otherwise mistake int's
    read-only ``real``/``imag``/``numerator``/``denominator`` descriptors (and the enum's own
    ``name``) for writable config fields and wrongly "enhance" the enum -- corrupting its repr and
    its pickling.  No genuine config struct subclasses int, so this is a safe, precise exclusion.
    """
    return isinstance(cls, type) and issubclass(cls, int) and cls not in (int, bool)


# ---------------------------------------------------------------------------
# enum pickling
# ---------------------------------------------------------------------------

def _reconstruct_enum(cls, value):
    """Pickle reconstructor for a Boost.Python enum value: rebuild the member from its int value."""
    return cls(value)


def _reduce_enum(member):
    return (_reconstruct_enum, (type(member), int(member)))


def _enable_enum_pickling(cls):
    """Make a Boost.Python enum picklable (its built-in ``__reduce_ex__`` is broken).

    Registered via ``copyreg`` so the pickler reconstructs the member from its int value; this also
    lets structs that *contain* enum fields (e.g. ``SolutionMetaData.endgame_success_code``) round-trip.
    """
    copyreg.pickle(cls, _reduce_enum)
    return cls


def _enhance_all(module):
    """Enhance every config-like class found in a bound module; make bound enums picklable."""
    for name in dir(module):
        obj = getattr(module, name)
        try:
            if _is_bound_enum(obj):
                _enable_enum_pickling(obj)
            elif _looks_like_config(obj):
                _enhance_config_class(obj)
        except Exception:
            # never let one odd member break importing the package
            pass


# ---------------------------------------------------------------------------
# owner (tracker / algorithm) helpers
# ---------------------------------------------------------------------------

def _config_holders(owner):
    """The objects whose configs count as this owner's, nearest first.

    A solver's settings are not all kept on the solver.  The tracker holds stepping, Newton and
    precision; the endgame holds the endgame, security and flavour configs.  Both are reached
    through live references (``get_tracker`` / ``get_endgame`` hand back the solver's own
    members), so writing through them writes into the solver -- which is what lets one flat
    ``solver.update(...)`` reach every setting the solve actually uses (b2#364).

    One level deep, deliberately: an endgame's own ``get_tracker()`` is the same tracker the
    solver would hand back, so recursing would add nothing but a second route to it.  Order is
    what breaks ties: the owner's own configs come first and win (see ``_field_owners``).
    """
    holders = [owner]
    for getter in ("get_tracker", "get_endgame"):
        fetch = getattr(owner, getter, None)
        if fetch is None:
            continue
        try:
            sub = fetch()
        except Exception:               # an owner that cannot produce one simply has none
            continue
        if sub is not None and _is_owner(sub):
            holders.append(sub)
    return holders


def _config_class_map(owner):
    """Map of {short_key: (holder, config_class)} over this owner and its config holders."""
    out = {}
    for holder in _config_holders(owner):
        for cls in holder.config_types():
            if cls is None:
                continue
            for key in (config_key(cls), cls.__name__.lower()):
                out.setdefault(key, (holder, cls))   # nearer holder wins
    return out


def _resolve_config_class(owner, key):
    """Find the (holder, config class) a config name refers to."""
    mapping = _config_class_map(owner)
    found = mapping.get(key) or mapping.get(str(key).lower())
    if found is None:
        names = config_names(owner)
        raise KeyError(
            "{0} has no config {1!r}.{2}  Available: {3}".format(
                type(owner).__name__, key, _did_you_mean(key, names), names))
    return found


def configure(self, **kwargs):
    """Change settings on this owner's configs in one call.

    Each keyword names a config (e.g. ``stepping``, ``newton``, ``tolerances``);
    its value is either a dict of fields to change, or a ready config object.

        tracker.configure(stepping={'max_step_size': 0.1},
                          newton={'max_num_newton_iterations': 2})

    Returns self.
    """
    for key, value in kwargs.items():
        holder, cls = _resolve_config_class(self, key)
        if isinstance(value, cls):
            cfg = value
        elif isinstance(value, dict):
            cfg = holder.get_config(cls).update(**value)
        else:
            raise TypeError(
                "value for {0!r} must be a dict of fields or a {1}, got {2}".format(
                    key, cls.__name__, type(value).__name__))
        holder.set_config(cfg)
    return self


def config_names(self):
    """The keyword names accepted by configure() for this owner, its tracker and its endgame."""
    return sorted({config_key(c)
                   for holder in _config_holders(self)
                   for c in holder.config_types() if c is not None})


def _field_owners(owner):
    """Map {field_name: (holder, config_class)} over this owner's configs, plus any collisions.

    Covers the configs of the owner and of the objects it keeps its settings in (see
    ``_config_holders``), so one flat call reaches a tracker's or an endgame's fields too.

    Two names for one field is resolved by distance, not refused: exactly one field name is
    shared across the library's config structs -- ``final_tolerance``, on both TolerancesConfig
    and EndgameConfig -- and on a solver the solver's own must keep winning, because
    ``solver.update(final_tolerance=...)`` is long-documented and reaches the endgame anyway
    through the solver's own push (b2#392).  A shadowed field stays reachable by naming its
    config: ``solver.configure(endgame={'final_tolerance': ...})``.

    A collision WITHIN one holder is still a collision, and update() refuses it rather than
    silently picking one.
    """
    owners = {}
    collisions = {}
    for holder in _config_holders(owner):
        local = {}
        for cls in holder.config_types():
            if cls is None:
                continue
            for f in writable_fields(cls):
                if f in local and local[f] is not cls:
                    collisions.setdefault(f, {local[f]})
                    collisions[f].add(cls)
                local[f] = cls
        for f, cls in local.items():
            owners.setdefault(f, (holder, cls))      # nearer holder wins
    return owners, collisions


def update(self, **fields):
    """Set config fields on this owner by NAME, each routed to whichever config owns it.

    You never name the config struct::

        solver.update(final_tolerance="1e-11",          # -> TolerancesConfig
                      max_num_crossed_path_resolve_attempts=3)   # -> ZeroDimConfig

    Strings work for every numeric field (converted exactly).  A field that none of this owner's
    configs has raises AttributeError with the valid names -- so a typo, or trying to set a tracker
    field on the algorithm (or vice versa), never silently does nothing.  Returns self, so calls
    chain.  To set a whole config at once, or to name the config explicitly, use configure().
    """
    owners, collisions = _field_owners(self)
    by_config = {}
    for key, value in fields.items():
        if key in collisions:
            raise AttributeError(
                "config field {0!r} is ambiguous on {1} (in {2}); set it with configure() naming "
                "the config".format(key, type(self).__name__,
                                    sorted(c.__name__ for c in collisions[key])))
        found = owners.get(key)
        if found is None:
            raise AttributeError(
                "{0} has no config field {1!r}.{2}  Valid fields: {3}".format(
                    type(self).__name__, key, _did_you_mean(key, owners), sorted(owners)))
        holder, cls = found
        by_config.setdefault((id(holder), cls), (holder, cls, {}))[2][key] = value
    # one get/update/set per touched config, not per field
    for holder, cls, kv in by_config.values():
        holder.set_config(holder.get_config(cls).update(**kv))
    return self


def get_settings(self, as_dict=False):
    """This owner's whole configuration as a carryable dict.

    By default (``as_dict=False``) the value is ``{config_name: config}`` -- each a copy of one of
    the owner's configs (e.g. ``{'stepping': SteppingConfig(...), 'tolerances': TolerancesConfig(...)}``),
    keyed by the short names config_names() lists.  The configs are independent copies (and picklable),
    so the dict is a plain Python value you can stash, tweak, and apply to other owners -- the way to
    carry one set of tracking settings across a series of related solves::

        settings = first_solver.get_settings()
        next_solver.set_settings(settings)

    Pass ``as_dict=True`` for the flat, human-readable view instead: one ``{field_name: value}`` dict
    across ALL configs (field names are unique across an owner's configs, so there is no collision).
    That form drops the struct layer nobody wants to poke at, and round-trips through ``set``::

        flat = solver.get_settings(as_dict=True)   # {'final_tolerance': ..., 'max_step_size': ..., ...}
        other.set(**flat)

    (The config-keyed default is kept because ``set_settings`` consumes it; use whichever fits.)  See
    set_settings() for applying the config-keyed form back.
    """
    configs = {}
    for holder in _config_holders(self):
        for cls in holder.config_types():
            if cls is not None:
                configs.setdefault(config_key(cls), holder.get_config(cls))
    if not as_dict:
        return configs
    flat = {}
    # nearest holder first, and the flat view must agree with update() about who owns a
    # shadowed name, so an earlier config is not overwritten by a later one
    for cfg in configs.values():
        for field, value in cfg.to_dict().items():
            flat.setdefault(field, value)
    return flat


def set_settings(self, settings, strict=False):
    """Apply a settings dict (from get_settings()) onto this owner.  Returns self.

    ``settings`` is ``{config_name: config}`` (or ``{config_name: {field: value}}``).  By default only
    the configs this owner actually has are applied and the rest are skipped -- so a bundle carried
    from one solver drops cleanly onto another whose config set differs (e.g. a different precision
    model, or a different algorithm stage).  Pass ``strict=True`` to instead raise on any key this
    owner does not have.
    """
    have = set(config_names(self))
    for key, value in dict(settings).items():
        if key not in have:
            if strict:
                raise KeyError(
                    "{0} has no config {1!r}.{2}  Available: {3}".format(
                        type(self).__name__, key, _did_you_mean(key, have), sorted(have)))
            continue
        holder, cls = _resolve_config_class(self, key)
        if isinstance(value, cls):
            holder.set_config(value)
        elif isinstance(value, dict):
            holder.set_config(holder.get_config(cls).update(**value))
        else:
            raise TypeError(
                "value for {0!r} must be a {1} or a dict of fields, got {2}".format(
                    key, cls.__name__, type(value).__name__))
    return self


def _accept_config_names(cls):
    """Let get_config()/set_config() take a config's short NAME as well as its class.

    ``config_names()`` advertises names such as ``'tolerances'``; before this, ``get_config`` took
    only the class and refused those very strings with a KeyError, so two adjacent methods on the
    same object did not compose (issue #392).  A name resolves through the same table
    ``configure()`` uses; the class form is unchanged.
    """
    native_get = cls.get_config
    native_set = cls.set_config

    def _holder_of(self, config_class):
        """Which of this object's config holders keeps the given config class, if any."""
        for holder in _config_holders(self):
            if any(c is config_class for c in holder.config_types()):
                return holder
        return None

    def get_config(self, config):
        """The named config, by class (``TolerancesConfig``) or by short name (``'tolerances'``,
        as ``config_names()`` lists them).

        Either may belong to this object's tracker or endgame rather than to the object itself
        (``solver.get_config('stepping')``, ``solver.get_config(EndgameConfig)``); it is fetched
        from whichever holds it, so every name ``config_names()`` lists can be asked for."""
        if isinstance(config, str):
            holder, config = _resolve_config_class(self, config)
        else:
            holder = _holder_of(self, config)
        if holder is not None and holder is not self:
            return holder.get_config(config)
        return native_get(self, config)

    def set_config(self, config):
        """Store a configuration struct, on this object or on whichever of its tracker and
        endgame owns that config's type."""
        holder = _holder_of(self, type(config))
        if holder is not None and holder is not self:
            return holder.set_config(config)
        return native_set(self, config)

    get_config.__doc__ = (native_get.__doc__ or "").rstrip() + "\n\n" + get_config.__doc__
    set_config.__doc__ = (native_set.__doc__ or "").rstrip() + "\n\n" + set_config.__doc__
    cls.get_config = get_config
    cls.set_config = set_config


def _enhance_owner_class(cls):
    """Attach configure()/config_names() to a tracker/algorithm class (idempotent)."""
    if getattr(cls, "_b2_owner_enhanced", False):
        return cls
    cls.configure = configure
    cls.config_names = config_names
    cls.update = update
    cls.set = update               # `set` reads more naturally to many users; same behavior
    cls.get_settings = get_settings
    cls.set_settings = set_settings
    if hasattr(cls, "get_config"):
        _accept_config_names(cls)
    cls._b2_owner_enhanced = True
    return cls


def _enhance_owners(module):
    """Attach owner helpers to every config-owning class in a bound module."""
    for name in dir(module):
        obj = getattr(module, name)
        try:
            if isinstance(obj, type) and _is_owner(obj):
                _enhance_owner_class(obj)
        except Exception:
            pass
