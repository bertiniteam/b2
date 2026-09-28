#!/usr/bin/env python
"""Extract and check the release notes for a tag: the top block of CHANGELOG.md.

The release notes are the block between the first two separator lines (79 underscores) of
CHANGELOG.md, HTML comments removed.  publish.yml calls this twice for a final release:

    python tools/release_notes.py --tag v4.0.0                       # check only
    python tools/release_notes.py --tag v4.0.0 --output TEMP_CHANGELOG.md

The first runs in the release's first job, before anything builds, so a CHANGELOG that is
not ready fails the release in seconds.  The second writes the notes the GitHub release
publishes.  A check fails when:

* the separator lines are missing, so there is no top block (the old one-line extraction
  returned an empty string here, and the release went out with an empty description);
* the block has no ``## [X.Y.Z] - YYYY-MM-DD`` heading, or its version is not the tag's;
* the heading's date is not a date -- typically ``unreleased``, left from development;
* the block has nothing under its heading;
* the block is longer than a GitHub release description may be.

Pull requests do not run this: whether the changelog names the next version correctly is a
release-time question.
"""

import argparse
import re
import sys
from pathlib import Path

SEPARATOR = "_" * 79                 # the line that fences each release's block
GITHUB_BODY_LIMIT = 125_000          # characters GitHub accepts in a release description
_HEADING = re.compile(r"^## \[(?P<version>[^\]]+)\] - (?P<date>.*?)\s*$", re.MULTILINE)
_DATE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


class ReleaseNotesError(Exception):
    """The changelog cannot supply release notes for the tag."""


def top_block(changelog_text):
    """The block between the first two separator lines, HTML comments removed.

    Parameters
    ----------
    changelog_text : str
        The whole of CHANGELOG.md.

    Returns
    -------
    str
        The block, starting at its opening separator line.
    """
    text = re.sub("<!--(.*?)-->", "", changelog_text, flags=re.DOTALL)
    start = text.find(SEPARATOR)
    if start < 0:
        raise ReleaseNotesError("CHANGELOG.md has no separator line (79 underscores), so there "
                                "is no top block to publish")
    end = text.find(SEPARATOR, start + len(SEPARATOR))
    if end < 0:
        raise ReleaseNotesError("CHANGELOG.md has only one separator line; the top block needs a "
                                "second one to end it")
    return text[start:end]


def check(block, tag):
    """Check that the block is publishable as the release notes for the tag.

    Parameters
    ----------
    block : str
        The top block, as returned by top_block.
    tag : str
        The release tag, e.g. ``v4.0.0`` (a leading ``v`` is optional).

    Returns
    -------
    tuple of (str, str)
        The heading's version and date.
    """
    heading = _HEADING.search(block)
    if heading is None:
        raise ReleaseNotesError("the top block has no '## [X.Y.Z] - YYYY-MM-DD' heading")
    version, date = heading.group("version"), heading.group("date")
    expected = tag[1:] if tag.startswith("v") else tag
    if version != expected:
        raise ReleaseNotesError(f"the top block is for {version}, but the tag is {tag}")
    if not _DATE.match(date):
        raise ReleaseNotesError(f"the top block's heading is dated {date!r}, not YYYY-MM-DD; "
                                f"date the {version} entry before tagging")
    if not block[heading.end():].strip():
        raise ReleaseNotesError(f"the {version} entry has nothing under its heading")
    if len(block) > GITHUB_BODY_LIMIT:
        raise ReleaseNotesError(f"the {version} entry is {len(block)} characters; a GitHub "
                                f"release description holds at most {GITHUB_BODY_LIMIT}")
    return version, date


def main(argv=None):
    """Check the release notes for a tag, and optionally write them out."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", required=True, help="the release tag, e.g. v4.0.0")
    parser.add_argument("--changelog", default="CHANGELOG.md", help="path to CHANGELOG.md")
    parser.add_argument("--output", help="write the release notes to this file")
    args = parser.parse_args(argv)

    try:
        block = top_block(Path(args.changelog).read_text())
        version, date = check(block, args.tag)
    except ReleaseNotesError as problem:
        print(f"::error title=Release notes::{problem}")
        return 1

    if args.output:
        Path(args.output).write_text(block)
    print(f"release notes for {args.tag}: {version}, dated {date}, {len(block)} characters")
    print(block)
    return 0


if __name__ == "__main__":
    sys.exit(main())
