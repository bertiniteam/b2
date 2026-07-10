"""Top-level ``bertini.randomize`` free function.

A thin, ergonomic wrapper over the :meth:`System.randomize` method so callers can write
``bertini.randomize(system, codimension)``.  The actual randomization logic lives in the C++ core.
"""


def randomize(system, codimension=None):
    """Randomize ``system``, returning a NEW system (the original is not mutated).

    When ``codimension`` is omitted the (overdetermined) system is squared down to its affine
    dimension.  Otherwise the system is randomized down to ``codimension`` functions -- a positive
    integer strictly less than the number of natural functions (``codimension=1`` yields a single
    function).
    """
    if codimension is None:
        return system.randomize()
    return system.randomize(codimension=codimension)
