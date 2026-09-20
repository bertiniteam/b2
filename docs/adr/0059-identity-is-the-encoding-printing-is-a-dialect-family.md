# 0059 — Identity is the canonical encoding; printing is a family of dialects, one per audience

**Status:** Accepted
**Date:** 2026-09-17

## Context

Two things were coupled that must not be.

**The content digest depended on the printer.**  Canonicalization orders the operands of a
sum or product by multidegree and breaks ties by comparing the operands' *printed* text
(`canonical.cpp`); folding like factors of a product matched bases by printed text too.
Consequences: (1) every cosmetic change to the printer moved the persistent content digest
(ADR-0042), so printing could not be improved without a `b2sysenc` bump; (2) the printer
walks the *tree*, so on a hash-consed graph a tiny DAG printed as an enormous expansion
(measured: 46 nodes, 196 KB), and building a large expression cost the expansion at every
level -- the second half of #417 (the first half, `MultiDegree`, was fixed by memoizing per
traversal, #418).

**One printer served three audiences.**  The single `print()` produced `^` for powers and
`(re,im)` for complex constants.  Bertini 1 reads `^` but *not* `(re,im)` -- its parser only
understands `(re+im*I)` -- so `to_classic_input()` wrote files Bertini 1 could not read, and
the Python layer grew a regex that rewrote `(re,im)` before parsing.  Python wants `**`, so
the bindings grew a character-replace shim for `repr`.  Neither shim knew the grammar it
was rewriting.

## Decision

1. **Identity is the canonical encoding, never printed text.**  The canonicalization
   tie-break and the like-factor key are the node's canonical encoding
   (`canonical_encoding.hpp`), the same exact text the digests are built on.  It is a pure
   function of content, deterministic, and linear in the DAG (a shared subtree is a
   back-reference after its first mention).  It is computed lazily -- only when two operands
   tie on multidegree -- and memoized for the sort.  The system encoding version is bumped
   to `b2sysenc/2`, because the encoder writes operands in canonicalization order and that
   order can differ where the old and new tie-breaks disagree.

2. **Printing is a family, one member per audience**, selected by a dialect carried on the
   output stream (`print_dialect.hpp`), so the recursive printers need no extra parameter
   and `operator<<` keeps working:
   - *Classic* (Bertini 1): `^`; a complex constant is `(re+im*I)` / `(re-im*I)`; every
     constant carries all of its digits.  This is `PrintClassic`, `System::to_classic_input`,
     the CLI's files, `node.to_classic()`, and the default of a stream nobody set a dialect on.
   - *PythonReadable* (`str`): `**`; otherwise the same constant shapes.
   - *PythonExact* (`repr`): `**`; constant spellings that `eval` in the `bertini` namespace
     rebuilds at full precision -- `real_mp('digits', precision)`, `Complex('re', 'im',
     precision)`, `Rational('p/q')`; integers and values a Python literal holds exactly stay
     bare.  `repr` is a fixed point under `eval` and the rebuilt values are bit-identical.
   Compatibility with another system is another dialect in this family, never a regex in
   Python.  The two shims are deleted.

## Consequences

- Do not reintroduce printed text into canonicalization, `Hash()`, or anything the digest
  depends on; `Node::Hash()` is per-process and printed text is presentation.  If a new
  ordering criterion is needed, derive it from the encoding.
- A printer change is now a presentation change: it never moves a digest and needs no
  `b2sysenc` bump.  The encoding, not the printer, is what must stay stable.
- Records written under `b2sysenc/1` are not recalled by `b2sysenc/2` builds (the version
  is part of every encoding and every ask).  This is the withdrawal of the "records read
  forever" promise noted in the 4.0.0 changelog; the reader (`System::FromCanonicalEncoding`)
  refuses a foreign version loudly rather than guessing.
- The complex spelling `(re,im)` is gone from every output.  Text produced by earlier
  versions that used it must be edited to `(re+im*I)` before Bertini 1 or this parser reads
  it; the parser never accepted the pair form itself.
- The reparse invariant of ADR-0012 is per dialect: Classic text reparses through the
  classic parser; PythonExact text rebuilds through `eval`.  Python-dialect text is not
  input for the classic parser (it does not read `**`).
