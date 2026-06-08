# ADR-0003: Linux wheel CI temporarily ran an import smoke test only (now reversed)

**Status:** Reversed 2026-06-08 — full Linux pytest restored and **confirmed green** on
x86_64 (run 27144247164: full suite `219 passed` inside the manylinux_2_34 container, all
platforms × py 3.10–3.14). The blocking crashes were fixed per ADR-0006 (uninitialized
slots / `dotfunc`) and ADR-0008 (writable-Ref corruption of `track_path`'s `end_time`).
Originally Accepted 2026-06-07.

> **Read this first (the short version).** For a period, the Linux wheel job ran only
> an `import bertini` smoke test instead of the full pytest suite, because the suite
> was crashing (SIGABRT/SIGSEGV) in CI. The crash was **first blamed on the container's
> old MPFR (3.1.6)**, but that was a **misdiagnosis**: the real cause is a
> version-independent bug — uninitialized `mpfr`/`mpc` numpy slots — now fixed in the
> bindings (**ADR-0006**). The crash reproduces on MPFR 4.x and even locally. The
> stopgap below has therefore been reversed; the full suite runs again.

## What actually happened (timeline)

1. **Original image:** `manylinux_2_28` (AlmaLinux 8, MPFR **3.1.6**). The full pytest
   suite SIGABRTed here. Because MPFR 3.1.6 has always-on `MPFR_ASSERTN`, the crash
   fired reliably, and it was *attributed* to MPFR 3.1.6.
2. **Stopgap (this ADR's original decision):** replace the Linux pytest run with an
   import smoke test so the build could go green, with the full suite still covered on
   macOS and Windows host runners.
3. **Image bumps that did NOT fix it:** the image was moved to `manylinux_2_34`
   (AlmaLinux 9, MPFR 4.1), then briefly `2_39`, and MPFR was even built from source —
   all in an attempt to "fix the SIGABRT" by getting a newer MPFR. None of it worked,
   because the bug is not MPFR-version dependent. Building MPFR from source was reverted
   (see commit `48afd397`); the image settled on **`manylinux_2_34`**, which is what CI
   uses today (`CIBW_MANYLINUX_X86_64_IMAGE: manylinux_2_34`).
4. **Real fix:** the uninitialized-slot bug was found and fixed in the eigenpy bindings
   — see **ADR-0006**. It is version-agnostic.
5. **Reversal (this ADR):** the full pytest suite was restored as
   `CIBW_TEST_COMMAND_LINUX` (with `CIBW_TEST_REQUIRES_LINUX: "pytest numpy"`), so it
   now runs inside the `manylinux_2_34` container — the in-container proof of the
   ADR-0006 fix.

> **Note on the merge artifact:** for a while `build_and_test.yml` and this ADR said
> "manylinux_2_28 / MPFR 3.1.6" while the live image was already `manylinux_2_34`. That
> contradiction came from a merge that combined the feature branch's `2_34` image with
> develop's smoke-test comments (written when the image was still `2_28`). Corrected
> 2026-06-08.

## Root cause (summary; full mechanism in ADR-0006)

An all-zero `mpfr_t`/`mpc_t` is Boost.Multiprecision's *uninitialized sentinel*, not a
valid zero. numpy zero-fills fresh user-dtype buffers (`NPY_NEEDS_INIT`), so any code
that **reads** a never-written slot and hands it to libmpfr/libmpc crashes (SIGSEGV);
and when the zero-fill guarantee is violated by malloc-dirty memory, a **write** onto a
garbage non-null `_mpfr_d` defeats BMP's null check and aborts (SIGABRT via
`MPFR_ASSERTN`). The MPFR version only affected *how reliably* the SIGABRT fired, not
whether the bug existed. ADR-0006 documents the three guards (getitem heal, setitem
zero-init, guarded ufunc/cast loops).

## Decision (reversed)

Run the full suite on Linux again:

```yaml
CIBW_TEST_REQUIRES_LINUX: "pytest numpy"
CIBW_TEST_COMMAND_LINUX: "cd {project} && python -m pytest python/test/ -q"
```

## Consequences

- **Full Linux coverage restored**, inside `manylinux_2_34` — the harshest realistic
  proof of the ADR-0006 fix (a real wheel, a real container, a real `import`).
- **Fallback if it regresses:** revert to the smoke test —
  `CIBW_TEST_COMMAND_LINUX: "python -c 'import bertini; print(bertini.__version__)'"`
  and drop `CIBW_TEST_REQUIRES_LINUX`. This is a last resort: a green smoke test would
  again hide real Python regressions on Linux.
- **Do not "fix" Linux test crashes by bumping MPFR or the manylinux image.** That was
  tried and does not address the root cause; the bindings guards (ADR-0006) do. Building
  MPFR from source in CI is specifically out of bounds.
