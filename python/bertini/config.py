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

#  individual authors of this file include:
#
#  silviana amethyst
#  UWEC
#

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

def _make_update(fields):
    def update(self, **kwargs):
        """Set one or more fields at once; returns self so calls can chain.

        Raises AttributeError on an unknown/misspelled field name.
        """
        for key, value in kwargs.items():
            if key not in fields:
                raise AttributeError(
                    "{0} has no config field {1!r}; valid fields: {2}".format(
                        type(self).__name__, key, list(fields)))
            setattr(self, key, value)
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


def enhance_config_class(cls):
    """Add update/to_dict/from_dict/repr/eq to a bound config class (idempotent)."""
    if getattr(cls, "_b2_config_enhanced", False):
        return cls
    fields = writable_fields(cls)
    cls._b2_fields = fields
    cls.update = _make_update(fields)
    cls.to_dict = _make_to_dict(fields)
    cls.from_dict = _make_from_dict()
    try:
        cls.__repr__ = _make_repr(fields)
        cls.__eq__ = _make_eq(fields)
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
            and len(writable_fields(cls)) > 0)


def enhance_all(module):
    """Enhance every config-like class found in a bound module."""
    for name in dir(module):
        obj = getattr(module, name)
        try:
            if _looks_like_config(obj):
                enhance_config_class(obj)
        except Exception:
            # never let one odd member break importing the package
            pass


# ---------------------------------------------------------------------------
# owner (tracker / algorithm) helpers
# ---------------------------------------------------------------------------

def _config_class_map(owner):
    """Map of {short_key: config_class} for the configs this owner accepts."""
    out = {}
    for cls in owner.config_types():
        if cls is None:
            continue
        out[config_key(cls)] = cls
        out[cls.__name__.lower()] = cls
    return out


def _resolve_config_class(owner, key):
    mapping = _config_class_map(owner)
    cls = mapping.get(key) or mapping.get(str(key).lower())
    if cls is None:
        valid = sorted({config_key(c) for c in owner.config_types() if c is not None})
        raise KeyError(
            "{0} has no config {1!r}; available: {2}".format(
                type(owner).__name__, key, valid))
    return cls


def configure(self, **kwargs):
    """Change settings on this owner's configs in one call.

    Each keyword names a config (e.g. ``stepping``, ``newton``, ``tolerances``);
    its value is either a dict of fields to change, or a ready config object.

        tracker.configure(stepping={'max_step_size': 0.1},
                          newton={'max_num_newton_iterations': 2})

    Returns self.
    """
    for key, value in kwargs.items():
        cls = _resolve_config_class(self, key)
        if isinstance(value, cls):
            cfg = value
        elif isinstance(value, dict):
            cfg = self.get_config(cls).update(**value)
        else:
            raise TypeError(
                "value for {0!r} must be a dict of fields or a {1}, got {2}".format(
                    key, cls.__name__, type(value).__name__))
        self.set_config(cfg)
    return self


def config_names(self):
    """The keyword names accepted by configure() for this owner."""
    return sorted({config_key(c) for c in self.config_types() if c is not None})


def enhance_owner_class(cls):
    """Attach configure()/config_names() to a tracker/algorithm class (idempotent)."""
    if getattr(cls, "_b2_owner_enhanced", False):
        return cls
    cls.configure = configure
    cls.config_names = config_names
    cls._b2_owner_enhanced = True
    return cls


def enhance_owners(module):
    """Attach owner helpers to every config-owning class in a bound module."""
    for name in dir(module):
        obj = getattr(module, name)
        try:
            if isinstance(obj, type) and _is_owner(obj):
                enhance_owner_class(obj)
        except Exception:
            pass
