# ADR-0066: A records session is the process

**Status:** Accepted
**Date:** 2026-09-27

## Context

Every recorded solve appends to a history file in `history/`, and each `OutputDirectory`
instance claims its own file on first append.  The name is
`YYYYMMDD_HHMMSS-pid<pid>[suffix].jsonl`, created exclusively ("one writer per file, ever",
`docs/records/b2rec-1.md`), and the suffix runs through 25 values: one process can claim at
most 25 names in any one second.

Solvers obtain their instance from `OutputDirectory::Shared(root)`, which already promised
"one instance -- hence one session history file -- per directory per process, however many
solvers attach", with a parameter sweep named as the reason.  But it held instances weakly.
A sweep's solvers attach, record, and release one after another, so between two solves
nothing held the instance; it was destroyed, and the next solve made a new one and claimed a
new file.  `parameter_sweep` over a grid of small systems starts more than 25 solves in a
second, and the 26th died with `OutputDirectory: could not claim a session history file`.
The parallel parameter homotopy tutorial's figure could not be regenerated because of it.

## Decision

`Shared()` holds instances strongly, for the life of the process.  A process records into a
directory through one instance and one history file, however many solves it runs.

The table is keyed by `weakly_canonical(absolute(root))`, and the instance keeps the absolute
path.  Without `absolute`, a relative path none of whose parts exist yet -- the default
`bertini_output` before its first solve creates it -- canonicalizes to itself, and after
creation to the absolute path, so a process's first recorded solve got a session of its own.
The stored absolute path also keeps a later `chdir` from moving an instance's writes.

Two cases replace the held instance with a fresh one: the directory has been deleted since
(a test clearing a fixed temporary path, a user clearing `bertini_output/`), and the caller is
a different process than the one that created it (a forked child must not write into its
parent's session file).  An instance whose session file is deleted while it lives claims a
new one rather than writing into the unlinked file.

Holding an instance must not mean holding its files.  A process may record into thousands of
directories -- the python test suite gives every test its own -- so on each `Shared()` call,
every instance only the table still holds closes its session file and results files, keeping
the session file's name and reopening that same file (append mode) on its next record; an
idle instance whose directory is gone is dropped from the table.  While in use, an instance
keeps at most `kMaxOpenResultsFiles` (16) per-run results files open, closing them all and
reopening on demand when the cap is reached.  Every stream is append-mode and flushed per
record, so closing one loses nothing.

## Consequences

- **Do not return `Shared()` to weak ownership** to "release" instances.  The instances are
  meant to outlive their solvers; a weak hold makes every solve in a sweep its own session,
  shreds the history across a file per solve, and crashes a fast sweep after 25 solves in a
  second.  `a_sweep_that_releases_between_solves_keeps_one_session` fails with the original
  error if the weak hold comes back.
- **Do not let a held instance keep its files open while idle, and do not let
  `results_streams_` grow without bound.**  Either runs a process out of file descriptors (256
  by default on macOS): one through many directories, the other through many runs in one.
  `many_directories_hold_no_idle_file_handles` counts the open handles on Linux (600 without
  the release, 0 with it), and `results_files_stay_correct_past_the_open_file_cap` covers the
  cap.
- **Do not key the table by `weakly_canonical(root)` alone.**  The default records path is
  relative and does not exist before the first solve.
  `a_relative_directory_is_one_session_from_its_first_solve` fails without `absolute`.
  A consequence: `records_path()` on a solver recording to the default directory returns the
  absolute path, not `bertini_output`.
- A process's history for a directory is one file, so a long session's file grows for as long
  as the process records.  Readers already scan every file in `history/` and do not care how
  records are divided among them.
- The fork and deleted-directory checks happen when a solver attaches through `Shared()`.  An
  instance obtained before a `fork()` and used directly in the child is not caught; attach
  through `Shared()` in the child.
