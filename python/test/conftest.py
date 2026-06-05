# This file is part of Bertini 2.
#
# python/test/conftest.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/test/conftest.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with python/test/conftest.py.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team
#
#  See <http://www.gnu.org/licenses/> for a copy of the license,
#  as well as COPYING.  Bertini2 is provided with permitted
#  additional terms in the b2/licenses/ directory.

#  individual authors of this file include:
#
#   silviana amethyst
#   spring 2025

"""Shared pytest fixtures for the Bertini Python test suite.

The multiprecision default precision is *global mutable state*
(``bertini.default_precision(n)``). Historically each ``unittest.TestCase`` set it in
``setUp``; a suite that forgot would inherit whatever an earlier-running suite left
behind, so a test could pass in a full run yet fail in isolation. The ``_reset_precision``
fixture below removes that hazard by restoring a known baseline before *every* test, so no
test can inherit a neighbor's precision.
"""

import pytest
import bertini as pb

DEFAULT_TEST_PRECISION = 30


@pytest.fixture(autouse=True)
def _reset_precision():
    """Reset the global default precision to a known baseline around every test.

    Autouse, so every test starts from ``DEFAULT_TEST_PRECISION`` regardless of run order,
    and the caller's precision is restored afterward.
    """
    old = pb.default_precision()
    pb.default_precision(DEFAULT_TEST_PRECISION)
    yield
    pb.default_precision(old)


@pytest.fixture
def precision(request):
    """Set the global default precision for a test, restoring it afterward.

    Defaults to ``DEFAULT_TEST_PRECISION``; parametrize indirectly to sweep precisions::

        @pytest.mark.parametrize("precision", [30, 50, 80], indirect=True)
        def test_something(precision):
            tol = mpfr_float(10) ** (-(precision - 3))
            ...
    """
    p = getattr(request, "param", DEFAULT_TEST_PRECISION)
    old = pb.default_precision()
    pb.default_precision(p)
    yield p
    pb.default_precision(old)
