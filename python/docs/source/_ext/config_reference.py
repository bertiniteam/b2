"""A ``config-reference`` directive: the settings reference, read out of the library.

The page that uses this directive is a few paragraphs of prose.  Every setting, default and
description on it comes from the code at build time, so adding a config field, changing a
default or rewriting a docstring updates the documentation with no documentation change --
which is the point: a hand-maintained table of defaults drifts from the code silently, and
that drift is what b2#406 was filed about.

What it finds, and how:

* **the config classes** -- every ``*Config`` class in the three modules that hold them.  The
  suffix is the library's own convention ("config classes live directly in the module, with
  Config-suffixed names"), so a new config appears here by being named like its siblings.
* **the fields** -- ``bertini.config.writable_fields``, the same introspection the settings
  surface itself routes by, so the page lists exactly what ``update()`` accepts.  It walks the
  MRO for data descriptors, seeing both plain members and property pairs, and skips read-only
  ones.
* **the defaults** -- by default-constructing the config and reading them.
* **the descriptions** -- each field's own docstring.

Two things it cannot read, and does not guess.  The adaptive-precision bounds are DERIVED from
a system rather than defaulted (see ``_DERIVED_FROM_SYSTEM``), and the tracker's settings that
have no config struct at all are hand-written on the page, because no amount of introspecting
config classes will find a method.
"""

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.statemachine import StringList
from sphinx.util.nodes import nested_parse_with_titles

#: Modules holding config classes, in the order the page presents them, with a line of prose
#: each.  A new module of configs is the one thing that needs adding here by hand.
_MODULES = [
    ("bertini.tracking", "Tracking", "How a path is followed: step sizes, Newton's method, and the precision model."),
    ("bertini.endgame", "Endgames", "How the last stretch to a singular endpoint is handled."),
    ("bertini.nag_algorithm", "Algorithms", "What the solver itself does, around the tracking."),
]

#: Fields whose default is not a value but a derivation.  A default-constructed ``AMPConfig``
#: leaves these unset -- they are meant to come from the system being tracked, through
#: ``bertini.tracking.amp_config_from(system)`` or ``set_bounds_from`` -- so printing whatever
#: an unset member happens to hold would be worse than saying where the value comes from.
_DERIVED_FROM_SYSTEM = {
    "coefficient_bound",
    "degree_bound",
    "linear_solve_error_bound",
    "jacobian_eval_error_bound",
    "function_eval_error_bound",
}


def _format_default(value):
    """Render a default value the way a reader would type it back in."""
    if isinstance(value, bool):                      # before int: a bool IS an int
        return "``{0}``".format(value)
    if type(value) is int or type(value) is float:
        return "``{0!r}``".format(value)
    if isinstance(value, int):
        # a bound enum -- Boost.Python makes these int subclasses, and str() gives either the
        # bare member name or the whole dotted path through the private module; neither is how
        # anyone writes it, so rebuild it as Type.Member
        return "``{0}.{1}``".format(type(value).__name__, str(value).rsplit(".", 1)[-1])
    text = str(value)
    return "``{0}``".format(text) if text else "(empty)"


def _describe(cls, field):
    """One line for a field, from its own docstring, with reST-hostile characters tamed."""
    doc = getattr(getattr(cls, field, None), "__doc__", None) or ""
    doc = " ".join(doc.split())
    # a bare identifier in a table cell can resolve to more than one target (the metadata
    # classes are exposed once per precision model); nothing here should become a link
    return doc.replace("|", r"\|") if doc else "*undocumented*"


def _config_classes(module):
    """The config classes a module exposes, by the library's own naming convention."""
    from bertini.config import writable_fields

    found = []
    for name in sorted(dir(module)):
        if not name.endswith("Config"):
            continue
        cls = getattr(module, name)
        if not isinstance(cls, type):
            continue
        if writable_fields(cls):
            found.append((name, cls))
    return found


def _rows_for(cls):
    """(field, default, description) for every settable field, or a reason there is none."""
    from bertini.config import writable_fields

    try:
        instance = cls()
    except Exception as e:                            # a config that cannot be default-built
        return [], "This config cannot be default-constructed here ({0}), so its defaults are not shown.".format(e)

    rows = []
    for field in sorted(writable_fields(cls)):
        if field in _DERIVED_FROM_SYSTEM:
            default = "*derived from the system*"
        else:
            try:
                default = _format_default(getattr(instance, field))
            except Exception:
                default = "*unavailable*"
        rows.append((field, default, _describe(cls, field)))
    return rows, None


class ConfigReference(Directive):
    """Emit the whole settings reference: every config class, field, default and description."""

    has_content = False

    def run(self):
        import importlib

        lines = []
        for module_name, heading, blurb in _MODULES:
            module = importlib.import_module(module_name)
            lines += [heading, "=" * len(heading), "", blurb, ""]

            for class_name, cls in _config_classes(module):
                qualified = "{0}.{1}".format(module_name, class_name)
                lines += [class_name, "-" * len(class_name), ""]
                lines += ["Set as ``{0}``, or by any of its field names directly.  "
                          "API: :class:`{1}`.".format(
                              _config_keyword(cls), qualified), ""]

                rows, problem = _rows_for(cls)
                if problem:
                    lines += [problem, ""]
                    continue

                lines += [".. list-table::", "   :header-rows: 1", "   :widths: 26 18 56", ""]
                lines += ["   * - Setting", "     - Default", "     - What it does"]
                for field, default, description in rows:
                    lines += ["   * - ``{0}``".format(field),
                              "     - {0}".format(default),
                              "     - {0}".format(description)]
                lines += [""]

        # nested_parse_with_titles, not plain nested_parse: the generated text has section
        # headings in it, and a directive's ordinary nested parse refuses those
        node = nodes.section()
        node.document = self.state.document
        nested_parse_with_titles(self.state, StringList(lines, source="<config-reference>"), node)
        return node.children


def _config_keyword(cls):
    """The short name ``configure()`` takes for this config."""
    from bertini.config import config_key

    return "{0}=...".format(config_key(cls))


def setup(app):
    app.add_directive("config-reference", ConfigReference)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
