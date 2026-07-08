"""ZeroDimConfig.recall (default True): an identical ask already in the records directory is recalled
rather than re-tracked.  recall=False forces a fresh track even when recorded -- the escape hatch for
path observers / benchmarking / re-verification.  (Bug context: a SolutionPathCollector silently
collected nothing on a re-solve, because the paths were recalled, not tracked.)
"""

import bertini as pb
from bertini import ZeroDimSolver
from bertini.nag_algorithm import ZeroDimConfig, observers as nobs


def _two_quadrics():
    x, y = pb.variables(['x', 'y'])
    sys = pb.System()
    sys.add_variable_group(x, y)
    sys.add_functions([x * x - 1, y * y - 1])     # 4 well-separated roots -> 4 total-degree paths
    return sys


def _paths_collected(rec_dir, recall):
    pb.random.set_random_seed(42)                 # SAME ask every call (system + settings + seed)
    solver = ZeroDimSolver(_two_quadrics())
    solver.record_to(str(rec_dir))                # isolate the records dir (no CWD pollution)
    cfg = solver.get_config(ZeroDimConfig)
    cfg.recall = recall
    solver.set_config(cfg)
    collector = nobs.SolutionPathCollector()
    solver.add_observer(collector)
    solver.solve()
    return len(collector.series)


def test_recall_default_true_then_false(tmp_path):
    rec = tmp_path / "records"
    assert _paths_collected(rec, recall=True) == 4    # 1st: fresh, tracks + records 4 paths
    assert _paths_collected(rec, recall=True) == 0    # 2nd: identical ask -> RECALLED, observer empty
    assert _paths_collected(rec, recall=False) == 4   # 3rd: recall=False -> forced fresh re-track


def test_recall_config_roundtrips():
    solver = ZeroDimSolver(_two_quadrics())
    assert solver.get_config(ZeroDimConfig).recall is True   # default
    cfg = solver.get_config(ZeroDimConfig)
    cfg.recall = False
    solver.set_config(cfg)
    assert solver.get_config(ZeroDimConfig).recall is False
