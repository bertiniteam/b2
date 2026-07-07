"""Tests for tools/assemble_versioned_docs.py.

Not collected by the package suite (pyproject testpaths = python/test); run explicitly:

    pytest tools/test_assemble_versioned_docs.py

These lock the accumulation behaviour so the multi-version docs deploy can be iterated on
locally in seconds, without a docs build.
"""

import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assemble_versioned_docs as avd  # noqa: E402


def make_site(tmp_path: Path, marker: str) -> Path:
    """A minimal built-site fixture with an underscore dir (to exercise .nojekyll)."""
    site = tmp_path / f"site_{marker}"
    (site / "python" / "_static").mkdir(parents=True, exist_ok=True)
    (site / "cpp").mkdir(exist_ok=True)
    (site / "python" / "index.html").write_text(f"<h1>{marker}</h1>")
    (site / "python" / "_static" / "t.css").write_text("x")
    (site / "cpp" / "index.html").write_text(f"<h1>{marker}</h1>")
    (site / "index.html").write_text("<html>landing</html>")
    (site / "style.css").write_text(":root{}")
    return site


def run(site, store, version, *, stable=False, date="2026-01-01"):
    argv = ["--site", str(site), "--store", str(store), "--version", version, "--date", date]
    if stable:
        argv.append("--stable")
    return avd.main(argv)


def load_json(store: Path):
    return json.loads((store / "versions.json").read_text())


def test_first_release_creates_everything(tmp_path):
    store = tmp_path / "store"
    rc = run(make_site(tmp_path, "3.0.0"), store, "3.0.0", stable=True, date="2026-07-14")
    assert rc == 0
    assert (store / "v3.0.0" / "python" / "index.html").read_text() == "<h1>3.0.0</h1>"
    assert (store / "stable" / "cpp" / "index.html").exists()
    assert (store / ".nojekyll").exists()
    assert (store / "style.css").exists()
    data = load_json(store)
    assert data["stable"] == "3.0.0"
    assert [v["version"] for v in data["versions"]] == ["3.0.0"]
    assert data["versions"][0]["stable"] is True
    assert "v3.0.0/" in (store / "index.html").read_text()


def test_newer_release_moves_stable_and_preserves_old(tmp_path):
    store = tmp_path / "store"
    run(make_site(tmp_path, "3.0.0"), store, "3.0.0", stable=True, date="2026-07-14")
    run(make_site(tmp_path, "3.1.0"), store, "3.1.0", stable=True, date="2026-09-01")

    # old version dir untouched, its recorded date preserved
    assert (store / "v3.0.0" / "python" / "index.html").read_text() == "<h1>3.0.0</h1>"
    data = load_json(store)
    assert data["stable"] == "3.1.0"
    assert [v["version"] for v in data["versions"]] == ["3.1.0", "3.0.0"]  # newest first
    dates = {v["version"]: v["released"] for v in data["versions"]}
    assert dates == {"3.0.0": "2026-07-14", "3.1.0": "2026-09-01"}
    # /stable/ now mirrors 3.1.0
    assert (store / "stable" / "python" / "index.html").read_text() == "<h1>3.1.0</h1>"


def test_patch_to_old_line_does_not_move_stable(tmp_path):
    store = tmp_path / "store"
    run(make_site(tmp_path, "3.0.0"), store, "3.0.0", stable=True, date="2026-07-14")
    run(make_site(tmp_path, "3.1.0"), store, "3.1.0", stable=True, date="2026-09-01")
    run(make_site(tmp_path, "3.0.1"), store, "3.0.1", stable=False, date="2026-09-15")

    data = load_json(store)
    assert data["stable"] == "3.1.0"  # unchanged
    assert [v["version"] for v in data["versions"]] == ["3.1.0", "3.0.1", "3.0.0"]
    assert (store / "stable" / "python" / "index.html").read_text() == "<h1>3.1.0</h1>"


def test_re_release_replaces_only_that_dir(tmp_path):
    store = tmp_path / "store"
    run(make_site(tmp_path, "3.0.0"), store, "3.0.0", stable=True)
    run(make_site(tmp_path, "3.1.0"), store, "3.1.0", stable=True)
    # rebuild 3.0.0 with new content; must not disturb 3.1.0
    run(make_site(tmp_path, "3.0.0-rebuilt"), store, "3.0.0", stable=False)
    assert (store / "v3.0.0" / "python" / "index.html").read_text() == "<h1>3.0.0-rebuilt</h1>"
    assert (store / "v3.1.0" / "python" / "index.html").read_text() == "<h1>3.1.0</h1>"


def test_versions_derived_from_dirs_present(tmp_path):
    """Truth is the directories; deleting one drops it from the derived views on next run."""
    store = tmp_path / "store"
    run(make_site(tmp_path, "3.0.0"), store, "3.0.0", stable=True)
    run(make_site(tmp_path, "3.1.0"), store, "3.1.0", stable=True)
    import shutil
    shutil.rmtree(store / "v3.0.0")
    # re-run for an existing version to regenerate the derived listing
    run(make_site(tmp_path, "3.1.0"), store, "3.1.0", stable=True)
    assert [v["version"] for v in load_json(store)["versions"]] == ["3.1.0"]


@pytest.mark.parametrize("bad", ["3.0.0rc1", "3.0.0.dev9", "3.0", "v3.0.0", "latest"])
def test_prerelease_or_malformed_version_rejected(tmp_path, bad):
    store = tmp_path / "store"
    assert run(make_site(tmp_path, "x"), store, bad) == 2
    assert not (store / "versions.json").exists()


def test_empty_site_rejected(tmp_path):
    store = tmp_path / "store"
    empty = tmp_path / "empty"
    empty.mkdir()
    assert run(empty, store, "3.0.0") == 2


def test_parse_version_ordering():
    keys = [avd.parse_version(n) for n in ("v3.0.0", "v3.1.0", "v3.0.1", "v10.0.0", "v3.0.0.post1")]
    assert None not in keys
    assert avd.parse_version("v3.1.0") > avd.parse_version("v3.0.9")
    assert avd.parse_version("v10.0.0") > avd.parse_version("v9.9.9")
    assert avd.parse_version("v3.0.0.post1") > avd.parse_version("v3.0.0")
    assert avd.parse_version("v3.0.0rc1") is None
