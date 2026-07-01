"""Bertini 2 tutorial: Evaluation of cyclic-n polynomials.

Assembles all code fragments from the tutorial into one runnable program.
Run:  python evaluation_cyclic.py
"""

import bertini
import numpy


def make_symbols(num_vars=10):
    """Programmatically create num_vars variables x0..x{n-1}."""
    x = [None] * num_vars  # preallocate the list
    for ii in range(num_vars):
        x[ii] = bertini.Variable('x' + str(ii))
    return x


def cyclic(vars):
    """Produce the cyclic-n polynomials from a list of variables."""
    n = len(vars)
    f = [None] * len(vars)
    y = []
    for ii in range(2):
        for x in vars:
            y.append(x)

    for ii in range(n - 1):
        f[ii] = numpy.sum([numpy.prod(y[jj:jj + ii + 1]) for jj in range(n)])

    # the last one is minus one
    f[-1] = numpy.prod(vars) - 1
    return f


def build_system(x):
    """Make a System and put the cyclic polynomials into it."""
    sys = bertini.System()

    for f in cyclic(x):
        sys.add_function(f)

    print(sys)  # long screen output, i know
    return sys


def associate_variables(sys, x):
    """Associate the variables with the system via a VariableGroup."""
    vg = bertini.VariableGroup()
    for var in x:
        vg.append(var)
    sys.add_variable_group(vg)


def simplify_system(sys):
    """Simplify the system in place (modifies the shared function tree)."""
    bertini.system.simplify(sys)


def evaluate_at_origin(sys):
    """Evaluate at the origin; all zeros except the last entry, which is -1."""
    s = numpy.zeros((10,), dtype=bertini.multiprec.Complex)
    result = sys.eval(s)
    assert complex(result[-1]) == -1                   # last cyclic function is (prod x) - 1
    assert all(complex(v) == 0 for v in result[:-1])   # the rest vanish at the origin
    return s


def evaluate_at_new_point(sys, s, num_vars=10):
    """Change the values of the vector and re-evaluate."""
    for ii in range(num_vars):
        s[ii] = bertini.multiprec.Complex(ii)
    result = sys.eval(s)
    return result


def main():
    num_vars = 10
    x = make_symbols(num_vars)
    sys = build_system(x)
    associate_variables(sys, x)
    simplify_system(sys)
    s = evaluate_at_origin(sys)
    evaluate_at_new_point(sys, s, num_vars)


if __name__ == '__main__':
    main()
