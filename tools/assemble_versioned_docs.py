#!/usr/bin/env python
"""Assemble a multi-version documentation tree for bertini2.org.

The docs site keeps *historical* versions: each real release lives forever under its own
``/vX.Y.Z/`` directory, a moving ``/stable/`` redirects to the newest release, and the site root is
a landing page that lists the versions.  This script performs the accumulation step: it takes a
freshly-built single-version ``site/`` and folds it into a persistent *store* directory (in CI, a
checkout of the ``docs-store`` branch), without ever rebuilding older versions.

GitHub Pages serves the store *directly from the ``docs-store`` branch* (Pages source = that
branch), so pushing the store IS the deploy -- there is no ``actions/deploy-pages`` step.  The
``/CNAME`` file (written when ``--cname`` is given) is what keeps the custom domain across builds.

Truth vs. derived (mirrors the records doctrine, ADR-0045/0047): the ``v*/`` directories present in
the store ARE the truth.  ``versions.json``, the root ``index.html``, and ``versions.html`` are
*derived, rebuildable views* -- regenerated from whatever version directories exist, every run.
Delete a ``v*/`` dir and re-run and it simply drops out of the listing.

The site root redirects straight to the current release -- most visitors do not want to pick a
version.  The human-readable chooser lives at ``/versions.html`` (root links to it), and
``/versions.json`` is the machine-readable index.

Layout produced in the store::

    /                     root: redirect to the current release (falls back to the chooser if none)
    /versions.html        the human-readable version chooser (lists every version)
    /versions.json        derived machine-readable version index
    /style.css            shared stylesheet (copied from the built site)
    /.nojekyll            so GitHub Pages serves _static/ etc. verbatim
    /CNAME                custom domain (only when --cname is given); persists it across builds
    /stable/              redirect stub to the newest release (a tiny page, NOT a byte-for-byte copy)
    /vX.Y.Z/              one durable directory per real release
    /vX.Y.Z/index.html      per-version landing (the built site's own index)
    /vX.Y.Z/{python,cpp,cli}/

Fast local loop (no CI, no compilation -- seconds)::

    mkdir -p /tmp/site/python /tmp/site/cpp
    echo hi > /tmp/site/python/index.html; echo hi > /tmp/site/cpp/index.html
    cp doc_resources/landing/index.html /tmp/site/index.html
    cp doc_resources/landing/style.css  /tmp/site/style.css
    python tools/assemble_versioned_docs.py --site /tmp/site --store /tmp/store --version 3.0.0 --stable
    python tools/assemble_versioned_docs.py --site /tmp/site --store /tmp/store --version 3.1.0 --stable
    python -m http.server -d /tmp/store 8000     # eyeball the whole site

Only real releases should be passed here; prereleases (.dev/a/b/rc) do not get a version directory.
"""

import argparse
import datetime
import html
import json
import re
import shutil
import sys
from pathlib import Path

# A version directory is 'v' + a PEP 440 public release, e.g. v3.0.0 (optionally v3.0.0.post1).
# Prereleases are intentionally NOT matched: they never get a durable directory.
_VDIR_RE = re.compile(r"^v(\d+)\.(\d+)\.(\d+)(?:\.post(\d+))?$")

STABLE_REDIRECT_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta http-equiv="refresh" content="0; url=/{vdir}/">
  <link rel="canonical" href="/{vdir}/">
  <title>Bertini 2 -- stable documentation</title>
</head>
<body>
  <p>The stable documentation is <a href="/{vdir}/">the latest release ({vdir})</a>.
     Redirecting&hellip;</p>
</body>
</html>
"""

# The site root sends visitors straight to the current release -- most people do not want to pick a
# version.  The human-readable version chooser lives at /versions.html (linked here for the few who
# do), and /versions.json is the machine-readable index.
ROOT_REDIRECT_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta http-equiv="refresh" content="0; url=/{vdir}/">
  <link rel="canonical" href="/{vdir}/">
  <title>Bertini 2 -- Documentation</title>
</head>
<body>
  <p>Redirecting to the <a href="/{vdir}/">latest documentation ({vdir})</a>&hellip;
     &nbsp;&middot;&nbsp; <a href="/versions.html">all versions</a></p>
</body>
</html>
"""

VERSIONS_PAGE_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Bertini 2 -- Documentation</title>
  <link rel="stylesheet" href="style.css">
  <style>
    .versions {{ list-style: none; margin: 0 0 3rem; padding: 0; }}
    .versions li {{
      display: flex; align-items: baseline; gap: 0.75rem;
      padding: 0.75rem 0; border-bottom: 1px solid var(--card-border);
    }}
    .versions li:first-child {{ border-top: 1px solid var(--card-border); }}
    .versions a {{ color: var(--accent); text-decoration: none; font-weight: 600; font-size: 1.1rem; }}
    .versions a:hover {{ text-decoration: underline; }}
    .versions .date {{ color: var(--muted); font-size: 0.9rem; margin-left: auto; }}
    .versions .tag {{
      font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.04em;
      color: var(--bg); background: var(--accent); border-radius: 999px; padding: 0.1rem 0.55rem;
    }}
  </style>
</head>
<body>
  <main>
    <header>
      <h1>Bertini 2</h1>
      <p class="subtitle">Numerical algebraic geometry &mdash; homotopy continuation for polynomial systems.</p>
    </header>

    <section class="cards">
      <a class="card" href="{stable_href}">
        <h2>Latest documentation</h2>
        <p class="lead">Stable release{stable_label}.</p>
        <p class="meta">Python, C++, and CLI documentation for the current release.</p>
        <span class="cta">Open the docs &rarr;</span>
      </a>
    </section>

    <h2>All versions</h2>
    <ul class="versions">
{version_items}
    </ul>

    <footer>
      <p>
        <a href="https://github.com/bertiniteam/b2">GitHub repository</a> &middot;
        <a href="https://pypi.org/project/bertini2/">PyPI package</a> &middot;
        Licensed under GPLv3 with additional terms.
      </p>
    </footer>
  </main>
</body>
</html>
"""


def parse_version(name: str):
    """Return a sort key tuple for a ``vX.Y.Z`` directory name, or ``None`` if it is not one."""
    m = _VDIR_RE.match(name)
    if not m:
        return None
    major, minor, patch, post = m.groups()
    return (int(major), int(minor), int(patch), int(post) if post else 0)


def copy_tree(src: Path, dst: Path) -> None:
    """Replace ``dst`` with a fresh copy of ``src`` (idempotent for re-releases).

    Sphinx build intermediates (``.doctrees``) are never published -- they are pure caches, so we
    drop them here as defense-in-depth even if the docs build forgot to redirect them out.
    """
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns(".doctrees"))


def write_stable_redirect(dst: Path, vdir: str) -> None:
    """Point ``/stable/`` at the newest release with a redirect stub (no byte-for-byte copy)."""
    if dst.exists():
        shutil.rmtree(dst)
    dst.mkdir(parents=True)
    (dst / "index.html").write_text(STABLE_REDIRECT_TEMPLATE.format(vdir=vdir))


def discover_versions(store: Path):
    """List the version directories present in ``store``, newest first."""
    found = []
    for child in store.iterdir():
        if child.is_dir():
            key = parse_version(child.name)
            if key is not None:
                found.append((child.name, key))
    found.sort(key=lambda t: t[1], reverse=True)
    return [name for name, _ in found]


def load_prior(store: Path):
    """Read the previous versions.json (dates + stable pointer) if present; else empty defaults."""
    path = store / "versions.json"
    if not path.exists():
        return {}, None
    try:
        data = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return {}, None
    dates = {v["version"]: v.get("released", "") for v in data.get("versions", [])}
    return dates, data.get("stable")


def render_versions_page(records, stable) -> str:
    """Build the /versions.html chooser HTML from the derived version records."""
    items = []
    for rec in records:
        tag = ' <span class="tag">stable</span>' if rec["version"] == stable else ""
        date = f'<span class="date">{html.escape(rec["released"])}</span>' if rec["released"] else ""
        items.append(
            f'      <li><a href="{html.escape(rec["path"])}">{html.escape(rec["version"])}</a>'
            f"{tag}{date}</li>"
        )
    stable_href = f"v{stable}/" if stable else (records[0]["path"] if records else "#")
    stable_label = f" (v{html.escape(stable)})" if stable else ""
    return VERSIONS_PAGE_TEMPLATE.format(
        stable_href=html.escape(stable_href),
        stable_label=stable_label,
        version_items="\n".join(items) if items else "      <li>No versions yet.</li>",
    )


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--site", type=Path, required=True,
                    help="freshly-built single-version site directory")
    ap.add_argument("--store", type=Path, required=True,
                    help="persistent accumulation directory (the docs-store checkout)")
    ap.add_argument("--version", required=True,
                    help="release version, e.g. 3.0.0 (directory becomes vX.Y.Z)")
    ap.add_argument("--stable", action="store_true",
                    help="also refresh /stable/ to mirror this version")
    ap.add_argument("--date", metavar="YYYY-MM-DD",
                    help="release date recorded for this version (default: today)")
    ap.add_argument("--style", type=Path,
                    help="stylesheet for the root landing (default: <site>/style.css if present)")
    ap.add_argument("--cname", metavar="DOMAIN",
                    help="write /CNAME with this custom domain (e.g. bertini2.org), so branch-source "
                         "GitHub Pages keeps the domain across builds")
    args = ap.parse_args(argv)

    site: Path = args.site
    store: Path = args.store

    if not site.is_dir() or not any(site.iterdir()):
        print(f"ERROR: --site is not a non-empty directory: {site}", file=sys.stderr)
        return 2
    if parse_version(f"v{args.version}") is None:
        print(f"ERROR: --version '{args.version}' is not a public release (X.Y.Z[.postN]); "
              "prereleases do not get a version directory.", file=sys.stderr)
        return 2

    store.mkdir(parents=True, exist_ok=True)
    (store / ".nojekyll").touch()
    if args.cname:
        (store / "CNAME").write_text(args.cname.strip() + "\n")

    vdir = f"v{args.version}"
    copy_tree(site, store / vdir)
    if args.stable:
        write_stable_redirect(store / "stable", vdir)

    # Shared stylesheet at the root, so the generated landing can reference /style.css.
    style_src = args.style or (site / "style.css")
    if style_src and style_src.is_file():
        shutil.copyfile(style_src, store / "style.css")

    # Derive the version list from the directories actually present (truth), carrying prior dates.
    prior_dates, prior_stable = load_prior(store)
    today = args.date or datetime.date.today().isoformat()
    records = []
    for name in discover_versions(store):
        version = name[1:]  # strip leading 'v'
        released = today if version == args.version else prior_dates.get(version, "")
        records.append({"version": version, "path": f"{name}/", "released": released})

    if args.stable:
        stable = args.version
    elif prior_stable and any(r["version"] == prior_stable for r in records):
        stable = prior_stable
    else:
        stable = records[0]["version"] if records else None
    for rec in records:
        rec["stable"] = rec["version"] == stable

    (store / "versions.json").write_text(
        json.dumps({"generated": today, "stable": stable, "versions": records}, indent=2) + "\n"
    )
    # Root sends visitors to the current release; the chooser is a separate /versions.html.  With no
    # stable version yet (bootstrap), fall back to serving the chooser at the root.
    (store / "versions.html").write_text(render_versions_page(records, stable))
    if stable:
        (store / "index.html").write_text(ROOT_REDIRECT_TEMPLATE.format(vdir=f"v{stable}"))
    else:
        (store / "index.html").write_text(render_versions_page(records, stable))

    print(f"OK: {vdir} written{' + stable redirect' if args.stable else ''}"
          f"{' + CNAME ' + args.cname if args.cname else ''}; "
          f"{len(records)} version(s) in store, stable=v{stable}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
