"""Tests for tools/py_doclint.py -- the wheel-free Python docstring linter.

Not collected by the package suite (pyproject testpaths = python/test); run explicitly:

    pytest tools/test_py_doclint.py

These lock the entity model, the numpydoc param-mismatch gate (pass 1), and the
undocumented ratchet (pass 2) so the linter can be iterated on in seconds without a build.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import py_doclint as pdl  # noqa: E402


def write_pkg(root: Path, files: dict) -> Path:
    """Materialise a fake package: {relative_path: source_text}. Returns the root."""
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)
    return root


def kinds(entities):
    return {e.qualname: e.kind for e in entities}


def by_name(entities):
    return {e.qualname: e for e in entities}


# --------------------------------------------------------------------------- #
# extract_documented_params
# --------------------------------------------------------------------------- #

def test_extract_params_numpydoc_basic():
    doc = (
        "Summary.\n\n"
        "Parameters\n----------\n"
        "functions : node or iterable\n"
        "    desc.\n"
        "variables : iterable of Variable\n"
        "    desc.\n"
    )
    assert pdl.extract_documented_params(doc) == {"functions", "variables"}


def test_extract_params_comma_and_star():
    doc = (
        "S.\n\nParameters\n----------\n"
        "x, y : int\n    d\n"
        "*args : tuple\n    d\n"
        "**kwargs : dict\n    d\n"
    )
    assert pdl.extract_documented_params(doc) == {"x", "y", "args", "kwargs"}


def test_extract_params_no_section_is_empty():
    assert pdl.extract_documented_params("Just a summary, no params.") == set()
    assert pdl.extract_documented_params(None) == set()


def test_extract_params_ignores_prose_colon():
    # Regression: an informal "Returns a solver: call .solve()" line (real case in
    # nag_algorithm.HomotopySolver) must NOT be read as a parameter named
    # "Returns a solver" -- names are identifiers, never multi-word prose.
    doc = (
        "S.\n\nParameters\n----------\n"
        "endgame : {'cauchy', 'powerseries'}\n\n"
        "Returns a solver: call ``.solve()`` then ``.all_solutions()``.\n"
    )
    assert pdl.extract_documented_params(doc) == {"endgame"}


def test_extract_params_stops_at_next_section():
    doc = (
        "S.\n\nParameters\n----------\n"
        "a : int\n    d\n"
        "Returns\n-------\n"
        "out : float\n    the result\n"
    )
    assert pdl.extract_documented_params(doc) == {"a"}


# --------------------------------------------------------------------------- #
# collect_entities
# --------------------------------------------------------------------------- #

def test_collect_public_private_dunder_and_nesting(tmp_path):
    src = (
        '"""Module doc."""\n'
        "def public_fn(a):\n"
        '    """doc"""\n'
        "    def inner():\n"          # nested -> ignored
        "        pass\n"
        "def _private_fn():\n"        # private -> ignored
        "    pass\n"
        "class Public:\n"
        '    """doc"""\n'
        "    def method(self, x):\n"
        '        """doc"""\n'
        "    def __init__(self):\n"   # dunder -> ignored
        "        pass\n"
        "    def _hidden(self):\n"    # private -> ignored
        "        pass\n"
        "class _Private:\n"          # private class -> ignored
        "    pass\n"
    )
    root = write_pkg(tmp_path / "pkg", {"m.py": src})
    ents = pdl.collect_entities(root)
    k = kinds(ents)
    assert k == {
        "m": "module",
        "m.public_fn": "function",
        "m.Public": "class",
        "m.Public.method": "method",
    }


def test_collect_init_module_name(tmp_path):
    root = write_pkg(tmp_path / "pkg", {"sub/__init__.py": '"""doc"""\ndef f():\n    """d"""\n'})
    names = {e.qualname for e in pdl.collect_entities(root)}
    assert "sub" in names          # not "sub.__init__"
    assert "sub.f" in names


def test_collect_skips_bin_and_pycache(tmp_path):
    files = {
        "real.py": '"""doc"""\n',
        "_bin/junk.py": "def x(): pass\n",
        "__pycache__/cached.py": "def y(): pass\n",
    }
    root = write_pkg(tmp_path / "pkg", files)
    names = {e.qualname for e in pdl.collect_entities(root)}
    assert names == {"real"}


def test_signature_params_strips_self(tmp_path):
    src = "class C:\n    def m(self, a, b, *c, **d):\n        pass\n"
    root = write_pkg(tmp_path / "pkg", {"m.py": src})
    m = by_name(pdl.collect_entities(root))["m.C.m"]
    assert m.params == ["a", "b", "c", "d"]


# --------------------------------------------------------------------------- #
# pass 1 -- param mismatch (regression)
# --------------------------------------------------------------------------- #

def test_param_mismatch_documented_absent_fails(tmp_path):
    src = (
        '"""m"""\n'
        "def f(a, b):\n"
        '    """S.\n\n    Parameters\n    ----------\n'
        "    a : int\n        d\n"
        "    zzz : int\n        not a real param\n"
        '    """\n'
    )
    root = write_pkg(tmp_path / "pkg", {"m.py": src})
    problems = pdl.find_param_mismatches(pdl.collect_entities(root))
    assert len(problems) == 1
    ent, bad = problems[0]
    assert ent.qualname == "m.f"
    assert bad == ["zzz"]


def test_param_mismatch_correct_docs_pass(tmp_path):
    src = (
        '"""m"""\n'
        "def f(a, b):\n"
        '    """S.\n\n    Parameters\n    ----------\n'
        "    a : int\n        d\n"
        "    b : int\n        d\n"
        '    """\n'
    )
    root = write_pkg(tmp_path / "pkg", {"m.py": src})
    assert pdl.find_param_mismatches(pdl.collect_entities(root)) == []


# --------------------------------------------------------------------------- #
# pass 2 -- undocumented ratchet + main() exit codes
# --------------------------------------------------------------------------- #

def _pkg_with_two_undocumented(tmp_path):
    src = (
        "def undoc_one():\n    pass\n"      # module has no docstring (undoc #1... plus func)
        "def undoc_two():\n    pass\n"
    )
    return write_pkg(tmp_path / "pkg", {"m.py": src})


def test_ratchet_equal_ok(tmp_path, capsys):
    root = _pkg_with_two_undocumented(tmp_path)
    n = len(pdl.find_undocumented(pdl.collect_entities(root)))
    baseline = tmp_path / "b.txt"
    baseline.write_text(str(n))
    rc = pdl.main(["--root", str(root), "--baseline", str(baseline)])
    assert rc == 0
    assert "at baseline" in capsys.readouterr().out


def test_ratchet_rose_fails(tmp_path):
    root = _pkg_with_two_undocumented(tmp_path)
    n = len(pdl.find_undocumented(pdl.collect_entities(root)))
    baseline = tmp_path / "b.txt"
    baseline.write_text(str(n - 1))
    rc = pdl.main(["--root", str(root), "--baseline", str(baseline)])
    assert rc == 1


def test_ratchet_dropped_ok(tmp_path):
    root = _pkg_with_two_undocumented(tmp_path)
    n = len(pdl.find_undocumented(pdl.collect_entities(root)))
    baseline = tmp_path / "b.txt"
    baseline.write_text(str(n + 5))
    rc = pdl.main(["--root", str(root), "--baseline", str(baseline)])
    assert rc == 0


def test_update_baseline_writes_current_count(tmp_path):
    root = _pkg_with_two_undocumented(tmp_path)
    n = len(pdl.find_undocumented(pdl.collect_entities(root)))
    baseline = tmp_path / "b.txt"
    rc = pdl.main(["--root", str(root), "--baseline", str(baseline), "--update-baseline"])
    assert rc == 0
    assert baseline.read_text().strip() == str(n)


def test_update_baseline_fails_on_broken_pass1(tmp_path):
    # A param mismatch must veto --update-baseline (never bless a baseline atop broken docs).
    src = (
        '"""m"""\n'
        "def f(a):\n"
        '    """S.\n\n    Parameters\n    ----------\n'
        "    ghost : int\n        d\n"
        '    """\n'
    )
    root = write_pkg(tmp_path / "pkg", {"m.py": src})
    baseline = tmp_path / "b.txt"
    rc = pdl.main(["--root", str(root), "--baseline", str(baseline), "--update-baseline"])
    assert rc == 1


def test_real_package_is_clean_at_committed_baseline():
    # The committed baseline must match the real package: pass 1 clean, pass 2 at baseline.
    rc = pdl.main([])
    assert rc == 0
