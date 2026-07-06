#!/usr/bin/env python
"""Release gate: fail if the tutorial timing numbers are stale.

The scaling-tutorial timing tables (``python/docs/source/tutorials/solving_at_scale/``) are
**re-measured on a reference machine**, not on CI -- wall-clock times are hardware- and
load-dependent, so CI can never regenerate them (see ``tools/update_scaling_timings.py`` and
``tools/refresh_doc_artifacts.py``).  That means CI cannot check them by regenerate-and-diff.  What
CI *can* do -- and what this script does -- is a cheap **provenance** check: the substitutions file
records the date the numbers were last measured (``.. |tw-date| replace:: YYYY-MM-DD``); if that is
older than a grace window, the release almost certainly went out without anyone re-running the
refresh on the reference box.

This is intended to run on the release-tag path in ``.github/workflows/publish.yml`` (not on every
PR).  On failure the fix is: on the reference machine, run

    python tools/refresh_doc_artifacts.py            # re-measures + rewrites the timing artifacts

commit the refreshed ``solving_at_scale/`` files, and re-tag.

Exit status: 0 fresh, 1 stale, 2 could not read/parse the date.
"""

import argparse
import datetime
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SUBS_FILE = REPO / "python" / "docs" / "source" / "tutorials" /  "parallelism" / "solving_at_scale" / "_timing_data.txt"

# Grace window.  Timings are refreshed at (roughly) each release; six months is comfortably longer
# than the release cadence, so tripping this means a refresh was almost certainly skipped, not that
# the numbers are merely a little old.  Tunable via --max-age-days.
DEFAULT_MAX_AGE_DAYS = 183

_TW_DATE = re.compile(r"^\.\.\s*\|tw-date\|\s*replace::\s*(\d{4}-\d{2}-\d{2})\s*$", re.MULTILINE)


def read_tw_date(subs_file: Path) -> datetime.date:
    """Extract the measured-on date from the substitutions file, or raise ValueError."""
    if not subs_file.exists():
        raise ValueError(f"timing substitutions file not found: {subs_file}")
    m = _TW_DATE.search(subs_file.read_text())
    if not m:
        raise ValueError(f"no `.. |tw-date| replace:: YYYY-MM-DD` line found in {subs_file}")
    return datetime.date.fromisoformat(m.group(1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--max-age-days", type=int, default=DEFAULT_MAX_AGE_DAYS,
                    help=f"fail if the timings are older than this many days (default {DEFAULT_MAX_AGE_DAYS})")
    ap.add_argument("--today", metavar="YYYY-MM-DD",
                    help="override the reference 'today' (for testing)")
    ap.add_argument("--subs-file", type=Path, default=SUBS_FILE,
                    help="path to the timing substitutions file (default: solving_at_scale/_timing_data.txt)")
    args = ap.parse_args(argv)

    try:
        measured = read_tw_date(args.subs_file)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 2

    today = datetime.date.fromisoformat(args.today) if args.today else datetime.date.today()
    age = (today - measured).days

    if age > args.max_age_days:
        print(f"STALE: tutorial timings were last measured {measured} ({age} days ago), "
              f"older than the {args.max_age_days}-day release grace window.", file=sys.stderr)
        print("Refresh them on the reference machine before releasing:", file=sys.stderr)
        print("    python tools/refresh_doc_artifacts.py", file=sys.stderr)
        print(f"then commit the updated {args.subs_file.parent.name}/ files and re-tag.", file=sys.stderr)
        return 1

    print(f"OK: tutorial timings measured {measured} ({age} days ago; "
          f"within the {args.max_age_days}-day window).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
