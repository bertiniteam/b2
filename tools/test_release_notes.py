"""Tests for tools/release_notes.py -- the release-time check on CHANGELOG.md's top block.

Not collected by the package suite (pyproject testpaths = python/test); run explicitly:

    pytest tools/test_release_notes.py
"""

import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import release_notes as rn  # noqa: E402

SEP = rn.SEPARATOR
REPO = Path(__file__).resolve().parent.parent


def changelog(top_heading, body="- a change\n", older="## [3.9.0] - 2026-01-02\n\n- older\n"):
    """A small changelog: a commented example, then two release blocks fenced by separators."""
    return (f"# Changelog\n\n<!--\n{SEP}\n\n## [0.0.0] - example\n-->\n\n"
            f"{SEP}\n\n{top_heading}\n\n{body}\n{SEP}\n\n{older}")


def test_extraction_matches_the_old_one_liner_on_the_real_changelog():
    # publish.yml used `text[start:text.find(sep, start+1)]` after stripping comments;
    # for a well-formed changelog the new extraction must publish exactly the same notes
    raw = (REPO / "CHANGELOG.md").read_text()
    text = re.sub("<!--(.*?)-->", "", raw, flags=re.DOTALL)
    start = text.find(SEP)
    old = text[start:text.find(SEP, start + 1)]
    assert rn.top_block(raw) == old


def test_a_dated_entry_for_the_tag_passes():
    block = rn.top_block(changelog("## [4.0.0] - 2026-10-01"))
    assert rn.check(block, "v4.0.0") == ("4.0.0", "2026-10-01")
    assert rn.check(block, "4.0.0") == ("4.0.0", "2026-10-01")   # the leading v is optional


def test_the_commented_example_is_never_the_top_block():
    block = rn.top_block(changelog("## [4.0.0] - 2026-10-01"))
    assert "0.0.0" not in block
    assert "4.0.0" in block and "3.9.0" not in block


def test_an_undated_entry_fails():
    block = rn.top_block(changelog("## [4.0.0] - unreleased"))
    with pytest.raises(rn.ReleaseNotesError, match="unreleased"):
        rn.check(block, "v4.0.0")


def test_a_version_other_than_the_tag_fails():
    block = rn.top_block(changelog("## [4.0.0] - 2026-10-01"))
    with pytest.raises(rn.ReleaseNotesError, match="tag is v4.0.1"):
        rn.check(block, "v4.0.1")


def test_a_missing_heading_fails():
    block = rn.top_block(changelog("# 4.0.0, no proper heading"))
    with pytest.raises(rn.ReleaseNotesError, match="heading"):
        rn.check(block, "v4.0.0")


def test_an_entry_with_nothing_under_its_heading_fails():
    block = rn.top_block(changelog("## [4.0.0] - 2026-10-01", body=""))
    with pytest.raises(rn.ReleaseNotesError, match="nothing under its heading"):
        rn.check(block, "v4.0.0")


def test_no_separator_fails_instead_of_publishing_empty_notes():
    with pytest.raises(rn.ReleaseNotesError, match="no separator"):
        rn.top_block("# Changelog\n\n## [4.0.0] - 2026-10-01\n\n- a change\n")


def test_one_separator_fails():
    with pytest.raises(rn.ReleaseNotesError, match="only one separator"):
        rn.top_block(f"# Changelog\n\n{SEP}\n\n## [4.0.0] - 2026-10-01\n\n- a change\n")


def test_a_separator_line_longer_than_79_still_fences_one_block():
    # the old extraction searched from start+1, so an 80-underscore line matched inside
    # itself and gave an empty block
    long_sep = SEP + "_"
    text = f"{long_sep}\n\n## [4.0.0] - 2026-10-01\n\n- a change\n{SEP}\n"
    assert "- a change" in rn.top_block(text)


def test_an_entry_longer_than_github_allows_fails():
    body = "- x\n" * (rn.GITHUB_BODY_LIMIT // 4 + 1)
    block = rn.top_block(changelog("## [4.0.0] - 2026-10-01", body=body))
    with pytest.raises(rn.ReleaseNotesError, match="at most"):
        rn.check(block, "v4.0.0")


def test_main_writes_the_notes_and_reports_errors_as_annotations(tmp_path, capsys):
    good = tmp_path / "CHANGELOG.md"
    good.write_text(changelog("## [4.0.0] - 2026-10-01"))
    out = tmp_path / "notes.md"
    assert rn.main(["--tag", "v4.0.0", "--changelog", str(good), "--output", str(out)]) == 0
    assert out.read_text() == rn.top_block(good.read_text())

    bad = tmp_path / "BAD.md"
    bad.write_text(changelog("## [4.0.0] - unreleased"))
    assert rn.main(["--tag", "v4.0.0", "--changelog", str(bad)]) == 1
    assert "::error title=Release notes::" in capsys.readouterr().out
