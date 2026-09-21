import contextlib
import os


def get_dll_paths():
    bertini_paths = os.getenv("BERTINI_WINDOWS_DLL_PATH")
    if bertini_paths is not None:
        return bertini_paths.split(os.pathsep)

    paths = []
    pkg_dir = os.path.dirname(__file__)

    # delvewheel bundles DLLs into bertini/.libs/
    paths.append(os.path.join(pkg_dir, ".libs"))

    # The package directory itself
    paths.append(pkg_dir)

    # Conda environment: CONDA_PREFIX/Library/bin
    conda_prefix = os.getenv("CONDA_PREFIX")
    if conda_prefix:
        paths.append(os.path.join(conda_prefix, "Library", "bin"))

    # Relative paths for when installed inside a conda environment
    # lib/python-version/site-packages/package -> ../../../../Library/bin
    paths.append(os.path.join(pkg_dir, "..", "..", "..", "..", "Library", "bin"))
    # Lib/site-packages/package -> ../../../Library/bin
    paths.append(os.path.join(pkg_dir, "..", "..", "..", "Library", "bin"))

    return paths


class DllDirectoryManager(contextlib.AbstractContextManager):
    """Restore DllDirectory state after importing Python module"""

    def add_dll_directory(self, dll_dir: str):
        # add_dll_directory can fail on relative path and non
        # existing path.
        # Since we don't know all the fail criterion we just ignore
        # thrown exception
        try:
            self.dll_dirs.append(os.add_dll_directory(dll_dir))
        except OSError:
            pass

    def __enter__(self):
        self.dll_dirs = []
        return self

    def __exit__(self, *exc_details):
        for d in self.dll_dirs:
            d.close()

    def __repr__(self):
        held = getattr(self, 'dll_dirs', None)
        if held is None:
            return '<DllDirectoryManager: not entered>'
        return '<DllDirectoryManager: %d directory%s added>' % (len(held),
                                                                '' if len(held) == 1 else 's')


def build_directory_manager():
    return DllDirectoryManager()
