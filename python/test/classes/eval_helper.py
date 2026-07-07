# This file is part of Bertini 2.  See b2/licenses/ for the license.

"""Test helper: evaluate a function-tree node at a point through the SLP.

Node-level evaluation was removed (the straight-line program is the sole evaluator), so tests that
used to set values on nodes and call node.eval_d()/eval_mp() now evaluate via node.eval(**point).
eval_at passes only the variables the node actually contains, so a single superset point serves a
function and all of its (variable-dropping) derivatives without tripping f.eval's typo guard.
"""


def eval_at(node, **point):
    names = {str(v) for v in node.variables()}
    return node.eval(**{name: value for name, value in point.items() if name in names})
