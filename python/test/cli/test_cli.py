# This file is part of Bertini 2.
#
# python/test/cli/test_cli.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/cli/test_cli.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/cli/test_cli.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

"""Packaging smoke test: the `bertini2` blackbox CLI is shipped in the wheel and runs.

This is an *interface/packaging* test, not a correctness test -- solver correctness is
gated by the C++ suite.  It proves that (1) the compiled CLI was installed into the
package, (2) the console-script shim resolves it, and (3) it loads the vendored shared
libs (--version reaches into libbertini2 for the version + dependency banner).

Skipped automatically when the bundled binary is absent (e.g. a from-source editable dev
checkout that never ran the wheel install step).
"""

import os
import shutil
import subprocess
import sys

import pytest

import bertini._cli as cli


def _cli_invocation():
    """How to invoke the CLI, or None if it is not available in this environment.

    Prefer the installed console-script (proves it landed on PATH); fall back to running
    the shim module directly when the bundled binary exists but PATH is not set up.
    """
    on_path = shutil.which("bertini2")
    if on_path:
        return [on_path]
    if os.path.exists(cli._binary_path()):
        return [sys.executable, "-m", "bertini._cli"]
    return None


pytestmark = pytest.mark.skipif(
    _cli_invocation() is None,
    reason="bundled bertini2 CLI not installed (from-source dev checkout)",
)


def _run(*args):
    return subprocess.run(
        _cli_invocation() + list(args),
        capture_output=True,
        text=True,
        timeout=60,
    )


def test_cli_version():
    result = _run("--version")
    assert result.returncode == 0, result.stderr
    assert "Bertini2" in result.stdout


def test_cli_help():
    result = _run("--help")
    assert result.returncode == 0, result.stderr
    assert "Usage" in result.stdout
