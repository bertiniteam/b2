#!/usr/bin/env python3
"""Check the classic emitter against a real Bertini 1.

Emits a small system with complex multiprecision, real multiprecision, rational and integer
coefficients through ``System.to_classic_input()``, runs Bertini 1 on the file in a scratch
directory, and matches every finite solution Bertini 1 reports against bertini2's own solve
of the same system.  Also checks that the emitted text reparses into bertini2 with the same
values.  This is the one check the test suite cannot carry: it needs a Bertini 1 install,
and the suites run with zero skips on every platform.

    python tools/bertini1_emitter_check.py                 # 'bertini' on PATH, temp scratch dir
    python tools/bertini1_emitter_check.py --bertini1 /path/to/bertini --scratch /tmp/b1run

Exit status 0 when every Bertini 1 solution matches a bertini2 solution and the counts
agree; 1 when Bertini 1 failed to run; 2 when the solutions disagree.  Bertini 1 writes its
many output files into the scratch directory, never into the working directory.
"""
import argparse
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

import bertini as pb
from bertini import Variable, System
from bertini.symbolics import Rational
from bertini.multiprec import complex_mp


def the_system():
    """Two quadrics in two variables, every constant kind the emitter spells: 4 solutions."""
    pb.default_precision(30)
    x, y = Variable('x'), Variable('y')
    c1 = pb.coefficient(complex_mp('0.3', '-0.5'))                       # complex, multiprecision
    c2 = pb.coefficient(complex_mp('1.2', '0.7'))
    third = pb.coefficient(complex_mp('0.333333333333333333333333333333', '0'))   # real, many digits
    s = System()
    s.add_variable_group([x, y])
    s.add_function(c1 * x**2 + y - Rational('6/5'))                     # a rational constant
    s.add_function(x + third * y**2 - c2 * x * y + 2)                   # an integer constant
    return s


def finite_solutions(path, num_vars):
    """Parse Bertini 1's finite_solutions: a count, then per solution a blank line and one
    'real imag' line per variable."""
    lines = path.read_text().split('\n')
    count = int(lines[0].strip())
    pairs = [tuple(float(t) for t in ln.split()) for ln in lines[1:] if ln.strip()]
    return [np.array([complex(*pairs[i * num_vars + k]) for k in range(num_vars)]) for i in range(count)]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--bertini1', default='bertini', help="the Bertini 1 binary (default: 'bertini' on PATH)")
    parser.add_argument('--scratch', help='directory Bertini 1 runs in (default: a fresh temporary directory)')
    parser.add_argument('--tolerance', type=float, default=1e-6, help='match tolerance, infinity norm (default 1e-6)')
    args = parser.parse_args(argv)

    scratch = pathlib.Path(args.scratch) if args.scratch else pathlib.Path(tempfile.mkdtemp(prefix='b2_bertini1_check_'))
    scratch.mkdir(parents=True, exist_ok=True)

    s = the_system()
    text = s.to_classic_input(mptype=2)
    (scratch / 'input').write_text(text)
    print(text[text.find('INPUT'):].rstrip())

    run = subprocess.run([args.bertini1, 'input'], cwd=scratch, capture_output=True, text=True, timeout=900)
    if run.returncode != 0:
        print(f'Bertini 1 exited with {run.returncode}:\n{run.stderr[-2000:]}')
        return 1

    b1_sols = finite_solutions(scratch / 'finite_solutions', s.num_variables())
    b2_sols = [np.array([complex(v) for v in sol]) for sol in pb.solve(s, seed=42, directory=str(scratch / 'b2_records'))]
    unmatched = [(sol, min(np.max(np.abs(sol - t)) for t in b2_sols)) for sol in b1_sols]
    unmatched = [(sol, d) for sol, d in unmatched if d > args.tolerance]

    back = pb.parse.system(text)
    pt = np.array([complex(0.37, -1.25), complex(2.5, 0.125)])
    residual = float(np.max(np.abs(np.asarray(back.eval(pt)) - np.asarray(s.eval(pt)))))

    print(f'Bertini 1: {len(b1_sols)} finite solutions; bertini2: {len(b2_sols)}; '
          f'{len(b1_sols) - len(unmatched)} of {len(b1_sols)} match to {args.tolerance:g}; '
          f'reparse residual {residual:g}; scratch {scratch}')
    for sol, d in unmatched:
        print(f'  unmatched Bertini 1 solution {sol}, nearest bertini2 solution at {d:g}')
    return 0 if not unmatched and len(b1_sols) == len(b2_sols) and residual == 0.0 else 2


if __name__ == '__main__':
    sys.exit(main())
