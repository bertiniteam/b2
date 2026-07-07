"""Console-script entry point for the bundled `bertini2` blackbox CLI.

The wheel ships the compiled solver at ``bertini/_bin/bertini2`` (see the
``install(TARGETS bertini2_exe ...)`` rule in ``core/CMakeLists.txt``).  The
``[project.scripts]`` entry in ``pyproject.toml`` maps the ``bertini2`` command to
``main`` below, which simply hands off to that binary -- so the launcher lands in the
environment's ``bin/`` (``Scripts\\bertini2.exe`` on Windows) uniformly across platforms.

This wheel CLI is threads-only.  MPI / multi-rank parallelism stays a from-source build.
"""

import os
import sys


def _binary_path():
    exe = "bertini2.exe" if os.name == "nt" else "bertini2"
    return os.path.join(os.path.dirname(__file__), "_bin", exe)


def main():
    exe = _binary_path()
    if not os.path.exists(exe):
        sys.stderr.write(
            "bertini2: bundled CLI not found at {}\n"
            "This wheel may have been built without the CLI.\n".format(exe)
        )
        return 1

    if os.name == "nt":
        # The exe needs the same vendored DLLs the extension module uses; prepend
        # their directories to PATH so the loader finds them.  Reuse the resolver
        # that windows_dll_manager already maintains for _pybertini.
        from .windows_dll_manager import get_dll_paths

        dll_paths = [p for p in get_dll_paths() if p]
        os.environ["PATH"] = os.pathsep.join(dll_paths + [os.environ.get("PATH", "")])
        # os.execv on Windows does not truly replace the process (the parent returns
        # immediately), which breaks shells waiting on the command.  Run as a child
        # and propagate the exit code instead.
        import subprocess

        return subprocess.run([exe, *sys.argv[1:]]).returncode

    # POSIX: replace this process with the solver so signals/exit codes pass through.
    os.execv(exe, [exe, *sys.argv[1:]])


if __name__ == "__main__":
    sys.exit(main())
