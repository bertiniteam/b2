"""bertini.random.complex_unit() must return a complex number of modulus exactly 1.

Regression: it used to normalize by sqrt(abs(z)) instead of abs(z), leaving modulus sqrt|z| (observed
magnitudes 0.76-1.09).  A unit-modulus complex is what the gamma trick needs (a well-scaled homotopy).
"""

import bertini


def test_complex_unit_has_modulus_one():
    for _ in range(25):
        z = complex(bertini.random.complex_unit())
        assert abs(abs(z) - 1.0) < 1e-12, abs(z)
