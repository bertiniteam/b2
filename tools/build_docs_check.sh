#!/usr/bin/env bash
# Local docs gate -- the SAME sphinx doctest + `-b html -W` (warnings-as-errors) build that CI
# (build_and_test's docs_gate job) and the release deploy (publish.yml) run.  Run it before pushing
# docstring / tutorial changes, so a reStructuredText mistake (e.g. a malformed inline literal) fails
# here rather than at the release deploy.
#
# Requires bertini importable in the current interpreter (`pip install -e .`, or a built wheel) --
# sphinx autodoc imports the package to render its docstrings.  The sphinx deps are installed below.
#
# Usage:  tools/build_docs_check.sh          # uses `python` on PATH
#         PYTHON=/path/to/python tools/build_docs_check.sh
set -euo pipefail

PY="${PYTHON:-python}"
cd "$(dirname "$0")/../python/docs"

echo "== installing sphinx deps =="
"$PY" -m pip install -q sphinx furo sphinxcontrib-bibtex gitpython matplotlib pandas sympy networkx

echo "== tutorial doctests (sphinx -b doctest) =="
"$PY" -m sphinx -b doctest --keep-going source ../../build/docs/doctest

echo "== html build, warnings as errors (sphinx -b html -W) =="
"$PY" -m sphinx -b html -W --keep-going -d ../../build/docs/py-doctrees source ../../build/docs/html

echo "OK: docs build clean -- doctests pass and no reST warnings."
