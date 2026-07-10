#!/usr/bin/env python3
#
# py_doclint.py -- lint the Bertini 2 pure-Python docstrings.
#
# The Python-side companion to tools/doclint.sh (which lints the C++ Doxygen
# docs).  Like that script it is deliberately cheap: it parses the pure-Python
# package with the standard-library `ast` module -- NO import of the compiled
# `_pybertini` extension, NO build, NO third-party linter -- so it finishes in a
# fraction of a second and can gate every PR in the same lane as the Doxygen job.
#
# Because it never imports anything, it sees only what the .py files *define*.
# The `from bertini._pybertini... import *` re-exports are opaque wildcards; the
# native classes/methods behind them live as Boost.Python docstrings in the C++
# `python_bindings/src/*_export.cpp` and are out of scope here (they need a build
# to observe).  This tool lints the ~28 files that actually define Python.
#
# It enforces two independent things, mirroring doclint.sh:
#
#   Pass 1 -- ERROR GATE (zero tolerance)
#     For every documented entity whose numpydoc docstring has a `Parameters`
#     section, each documented parameter name must exist in the actual signature.
#     A documented-but-absent name means the docs describe code that does not
#     exist (the analogue of Doxygen's @param-name-mismatch warning).  MISSING
#     param docs are NOT flagged here (matches Doxygen WARN_NO_PARAMDOC = NO), so
#     the undocumented backlog never blocks a build.  Any hit fails the lint.
#
#   Pass 2 -- UNDOCUMENTED RATCHET (monotonic)
#     Count public entities (module, public class, public function, public
#     method) that lack a docstring, and compare against
#     tools/py_doc_undocumented_baseline.txt.  The count may only decrease: if it
#     rises the lint fails; when it drops, run with --update-baseline to lock in
#     the gain.  Target is 0, at which point the gate can be made absolute.
#
# "Public" = a name that does not start with an underscore (dunders such as
# __init__ are therefore NOT required to carry their own docstring -- the class
# docstring covers them).  Only module-level defs and one level of methods are
# considered; nested/inner defs are ignored.  A file named with a leading
# underscore (e.g. _calculus.py) is still linted -- it defines public functions
# that the package re-exports.
#
# Usage:
#   tools/py_doclint.py                  # run both passes (CI mode)
#   tools/py_doclint.py --update-baseline  # rewrite the baseline from current count
#   tools/py_doclint.py --root <dir>     # lint a different package root (tests)
#
# Exit status: 0 if both passes pass (or --update-baseline succeeds atop a clean
# pass 1), non-zero otherwise.

import argparse
import ast
import os
import re
import sys

# Package root and baseline default to this repo's layout; overridable for tests.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ROOT = os.path.join(_REPO_ROOT, "python", "bertini")
DEFAULT_BASELINE = os.path.join(_REPO_ROOT, "tools", "py_doc_undocumented_baseline.txt")


class Entity:
    """One documentable Python entity: a module, class, function, or method."""

    def __init__(self, kind, qualname, path, lineno, has_doc, params, docstring):
        self.kind = kind            # 'module' | 'class' | 'function' | 'method'
        self.qualname = qualname    # e.g. "records.OutputDirectory.rebuild"
        self.path = path            # repo-relative source path
        self.lineno = lineno        # 1-based line of the def (0 for a module)
        self.has_doc = has_doc      # bool: a non-empty docstring is present
        self.params = params        # list[str]: actual signature parameter names
        self.docstring = docstring  # str | None: the raw docstring, if any

    def __repr__(self):
        return "Entity(%s %s doc=%s)" % (self.kind, self.qualname, self.has_doc)


def _is_public(name):
    """True for a name that does not start with an underscore (dunders excluded)."""
    return not name.startswith("_")


def _signature_params(node):
    """The parameter names of a def node, stripped of self/cls and */** markers."""
    a = node.args
    names = []
    for arg in list(a.posonlyargs) + list(a.args) + list(a.kwonlyargs):
        names.append(arg.arg)
    if a.vararg:
        names.append(a.vararg.arg)
    if a.kwarg:
        names.append(a.kwarg.arg)
    # self/cls are never documented as parameters; drop a leading one.
    if names and names[0] in ("self", "cls"):
        names = names[1:]
    return names


def _docstring_of(node):
    """(has_doc, raw_docstring) for a module/class/def node."""
    doc = ast.get_docstring(node, clean=False)
    if doc is not None and doc.strip():
        return True, doc
    return False, None


_SECTION_HEADERS = {
    "parameters", "returns", "yields", "raises", "examples", "notes",
    "see also", "references", "attributes", "other parameters", "warnings",
    "warns", "methods",
}


def extract_documented_params(docstring):
    """Parameter names declared in a numpydoc ``Parameters`` section.

    Conservative on purpose (Pass 1 is zero-tolerance and cannot be ratcheted):
    only ``name : type`` entry lines at the section's base indentation are read,
    so prose and continuation lines never masquerade as parameters.  Returns the
    set of declared names (``*args``/``**kwargs`` markers stripped).  A docstring
    with no ``Parameters`` section yields an empty set.
    """
    if not docstring:
        return set()
    lines = docstring.splitlines()
    # Find "Parameters" followed by an underline of dashes.
    start = None
    for i in range(len(lines) - 1):
        if lines[i].strip().lower() == "parameters" and set(lines[i + 1].strip()) == {"-"}:
            start = i + 2
            break
    if start is None:
        return set()

    names = set()
    base_indent = None
    # A numpydoc entry line is `name : type` or `name1, name2 : type`.  Each name
    # is a (possibly *-prefixed) identifier -- NO internal spaces -- so informal
    # prose like "Returns a solver: call .solve()" never matches as a parameter.
    entry_re = re.compile(
        r"^(?P<indent>[ \t]*)"
        r"(?P<names>\*{0,2}\w+(?:\s*,\s*\*{0,2}\w+)*)"
        r"\s*:\s+\S")
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped:
            continue
        indent = len(line) - len(line.lstrip())
        # A new section ("Returns" etc.) ends the Parameters block.
        if stripped.lower() in _SECTION_HEADERS:
            break
        if base_indent is None:
            base_indent = indent
        # Only lines at the base indent are parameter entries; deeper lines are
        # the entry descriptions.  A dedent below the base indent ends the block.
        if indent < base_indent:
            break
        if indent > base_indent:
            continue
        m = entry_re.match(line)
        if not m:
            continue
        for raw in m.group("names").split(","):
            name = raw.strip().lstrip("*").strip()
            if name.isidentifier():
                names.add(name)
    return names


def collect_entities(root):
    """All documentable entities under ``root``, in deterministic (path, line) order.

    Skips ``_bin/`` and ``__pycache__/``.  For each ``.py`` file yields the module
    entity plus its public top-level classes/functions and one level of public
    methods.  Nested and inner defs are ignored.
    """
    entities = []
    py_files = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in ("_bin", "__pycache__")]
        for fn in filenames:
            if fn.endswith(".py"):
                py_files.append(os.path.join(dirpath, fn))
    py_files.sort()

    for path in py_files:
        rel = os.path.relpath(path, root)
        with open(path, "r", encoding="utf-8") as fh:
            source = fh.read()
        tree = ast.parse(source, filename=path)

        has_doc, doc = _docstring_of(tree)
        modname = rel[:-3].replace(os.sep, ".")
        # `pkg/__init__.py` names the package itself, not a `pkg.__init__` module.
        if modname.endswith(".__init__"):
            modname = modname[: -len(".__init__")]
        elif modname == "__init__":
            modname = "(package root)"
        entities.append(Entity("module", modname, rel, 0, has_doc, [], doc))

        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and _is_public(node.name):
                hd, d = _docstring_of(node)
                entities.append(Entity(
                    "function", "%s.%s" % (modname, node.name), rel,
                    node.lineno, hd, _signature_params(node), d))
            elif isinstance(node, ast.ClassDef) and _is_public(node.name):
                hd, d = _docstring_of(node)
                entities.append(Entity(
                    "class", "%s.%s" % (modname, node.name), rel,
                    node.lineno, hd, [], d))
                for sub in node.body:
                    if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)) and _is_public(sub.name):
                        mhd, md = _docstring_of(sub)
                        entities.append(Entity(
                            "method", "%s.%s.%s" % (modname, node.name, sub.name), rel,
                            sub.lineno, mhd, _signature_params(sub), md))
    return entities


def find_param_mismatches(entities):
    """(entity, sorted[bad_names]) for every entity documenting an absent parameter."""
    problems = []
    for e in entities:
        if e.kind not in ("function", "method") or not e.has_doc:
            continue
        documented = extract_documented_params(e.docstring)
        if not documented:
            continue
        bad = documented - set(e.params)
        if bad:
            problems.append((e, sorted(bad)))
    return problems


def find_undocumented(entities):
    """Entities lacking a docstring, in collection order."""
    return [e for e in entities if not e.has_doc]


def _read_baseline(path):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            text = fh.read().strip()
        return int(text) if text else None
    except (OSError, ValueError):
        return None


def main(argv=None):
    parser = argparse.ArgumentParser(description="Lint pure-Python docstrings (wheel-free).")
    parser.add_argument("--root", default=DEFAULT_ROOT,
                        help="package root to lint (default: python/bertini)")
    parser.add_argument("--baseline", default=DEFAULT_BASELINE,
                        help="undocumented-count baseline file")
    parser.add_argument("--update-baseline", action="store_true",
                        help="rewrite the baseline from the current count")
    args = parser.parse_args(argv)

    entities = collect_entities(args.root)

    # ---- Pass 1: correctness (documented params must exist) ----------------
    problems = find_param_mismatches(entities)
    print("== py doclint pass 1: documentation correctness ==")
    if problems:
        print("FAIL: %d entit(ies) document a parameter that does not exist:"
              % len(problems), file=sys.stderr)
        print(file=sys.stderr)
        for e, bad in problems:
            print("  %s:%d: %s documents unknown param(s): %s"
                  % (e.path, e.lineno, e.qualname, ", ".join(bad)), file=sys.stderr)
        pass1_ok = False
    else:
        print("OK: no documentation correctness errors.")
        pass1_ok = True

    # ---- Pass 2: undocumented ratchet --------------------------------------
    undoc = find_undocumented(entities)
    undoc_count = len(undoc)

    if args.update_baseline:
        with open(args.baseline, "w", encoding="utf-8") as fh:
            fh.write("%d\n" % undoc_count)
        print("== baseline updated: %d undocumented entit(ies) ==" % undoc_count)
        # Never bless a baseline atop broken docs.
        return 0 if pass1_ok else 1

    baseline = _read_baseline(args.baseline)
    print("== py doclint pass 2: undocumented ratchet ==")
    if baseline is None:
        print("note: no baseline file; current undocumented count is %d." % undoc_count)
        print("      run 'tools/py_doclint.py --update-baseline' to seed it.")
        pass2_ok = True
    elif undoc_count > baseline:
        print("FAIL: undocumented entities rose from %d to %d (+%d)."
              % (baseline, undoc_count, undoc_count - baseline), file=sys.stderr)
        print("      document the new code, or this is genuinely new public surface:",
              file=sys.stderr)
        for e in undoc:
            print("  %s:%d: %s %s" % (e.path, e.lineno, e.kind, e.qualname), file=sys.stderr)
        pass2_ok = False
    else:
        if undoc_count < baseline:
            print("OK: undocumented entities dropped from %d to %d (-%d)."
                  % (baseline, undoc_count, baseline - undoc_count))
            print("    run 'tools/py_doclint.py --update-baseline' to lock in the gain.")
        else:
            print("OK: undocumented entities at baseline (%d)." % undoc_count)
        pass2_ok = True

    return 0 if (pass1_ok and pass2_ok) else 1


if __name__ == "__main__":
    sys.exit(main())
