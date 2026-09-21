# This file is part of Bertini 2.
#
# python/bertini/_repr.py is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# python/bertini/_repr.py is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this file.  If not, see <http://www.gnu.org/licenses/>.
#
#  Copyright(C) Bertini2 Development Team

"""Make printing a bertini object tell you something about it.

Without this, most of the library prints the default
``<bertini._pybertini.tracking.AMPTracker object at 0x7f3c...>`` -- the type you already
knew and an address you cannot use.  This module replaces that for every public type that
had it, by attaching ``__repr__`` at import; ``str()`` follows ``repr()`` for any type
that does not define its own, so one function serves both.

Two shapes, chosen by what the object is:

* **Machinery** -- trackers, endgames, observers, collectors, decompositions -- reads
  ``<AMPTracker: homotopy in 2 variables, RKF45, tolerance 1e-05, adaptive precision at
  30 digits, at t=0.31 after 47 steps>``.  The angle brackets are the honest signal that
  the text is a description and not a way to rebuild the object; none of these could be
  reconstructed from any amount of text, and most hold a C++ object that a python
  expression cannot name.
* **Values** -- configs, metadata, ``StraightLineHomotopy`` -- already read
  ``SteppingConfig(initial_step_size=..., ...)`` and keep that shape, which is closer to
  something you could type back in.

A tracker's description is **live**: it reports where the path is now, so printing one
inside an observer callback -- the moment you most want to look -- says where tracking
has got to rather than repeating its settings.

Nothing here may raise.  A repr that throws turns an ordinary ``print`` in a debugging
session into a traceback about printing, and every accessor used below reaches into C++
state that may not be there yet, so each one goes through :func:`_ask`.
"""

def _ask(getter, default=None):
    """Call an accessor, yielding `default` if it is absent or unhappy."""
    try:
        value = getter()
    except Exception:
        return default
    return default if value is None else value


def _number(value, default='?'):
    """A short decimal for a number of any of the multiprecision or builtin types."""
    try:
        as_float = float(value)
    except Exception:
        return default
    if as_float == int(as_float) and abs(as_float) < 1e16:
        return '%d' % int(as_float)
    return '%g' % as_float


def _time(value, default='?'):
    """A time, as `0.31` when it is real and `0.31+0.5i` when it is not."""
    try:
        real, imaginary = float(value.real), float(value.imag)
    except Exception:
        return default
    if imaginary == 0.0:
        return _number(real)
    return '%s%s%si' % (_number(real), '+' if imaginary > 0 else '-', _number(abs(imaginary)))


def _count(number, noun, default=None):
    """`1 variable`, `3 variables`, or `default` when the number is unavailable."""
    if number is None:
        return default
    return '%d %s%s' % (number, noun, '' if number == 1 else 's')


def _described(obj, description):
    """The angle-bracket form: the type, and what there is to say about it."""
    name = type(obj).__name__
    return '<%s: %s>' % (name, description) if description else '<%s>' % name


def _has_its_own_repr(cls):
    """Does this type get a __repr__ from anywhere but the default, or from this module?

    A description installed here does not count as the type's own: a collector inherits
    from an observer, so the bare name attached to the base must not stand in the way of
    the count attached to the collector itself.
    """
    for ancestor in cls.__mro__:
        if '__repr__' in ancestor.__dict__:
            if ancestor.__name__ in ('object', 'instance'):
                return False
            return getattr(ancestor.__dict__['__repr__'], '__module__', None) != __name__
    return False


def _attach(cls, function, replacing=False):
    """Give `cls` this __repr__.

    Only fills in for the default one, so a description written anywhere else -- in the
    bindings, or in the python class itself -- wins over this module.  `replacing` is for
    the few types whose existing repr is the thing being fixed.
    """
    if cls is None:
        return
    if not replacing and _has_its_own_repr(cls):
        return
    try:
        cls.__repr__ = function
    except (TypeError, AttributeError):
        pass       # a type that refuses the assignment keeps the default; not worth failing over


def _each(classes, function, replacing=False):
    for cls in classes:
        _attach(cls, function, replacing)


def _classes(module, names):
    """The named classes of a module, skipping any that are not there."""
    return [getattr(module, name) for name in names if hasattr(module, name)]


# ----- trackers ---------------------------------------------------------------------

_PRECISION_MODE = {
    'AMPTracker': 'adaptive precision',
    'DoublePrecisionTracker': 'double precision',
    'MultiplePrecisionTracker': 'fixed precision',
}


def _tracker_repr(self):
    parts = []

    system = _ask(self.get_system)
    if system is not None:
        kind = 'homotopy' if _ask(system.have_path_variable, False) else 'system'
        described = _count(_ask(system.num_variables), 'variable')
        parts.append('%s in %s' % (kind, described) if described else kind)

    predictor = _ask(self.predictor)
    if predictor is not None:
        parts.append(str(predictor).rsplit('.', 1)[-1])

    tolerance = _ask(self.tracking_tolerance)
    if tolerance is not None:
        parts.append('tolerance %s' % _number(tolerance))

    mode = _PRECISION_MODE.get(type(self).__name__, 'precision')
    digits = _ask(self.current_precision)
    parts.append('%s at %d digits' % (mode, digits) if digits else mode)

    point = _ask(self.current_point)
    if point is None or len(point) == 0:
        parts.append('not started')
    else:
        steps = _ask(self.num_total_steps_taken)
        where = 'at t=%s' % _time(_ask(self.current_time))
        parts.append('%s after %s' % (where, _count(steps, 'step')) if steps else where)

    return _described(self, ', '.join(parts))


# ----- endgames ---------------------------------------------------------------------

def _endgame_repr(self):
    parts = []

    boundary = _ask(self.boundary_time)
    if boundary is not None:
        parts.append('from t=%s' % _time(boundary))

    target = _ask(self.target_time)
    if target is not None:
        parts.append('to t=%s' % _time(target))

    cycle = _ask(self.cycle_number)
    approximation = _ask(self.final_approximation)
    if approximation is None or len(approximation) == 0:
        parts.append('not yet run')
    elif cycle:
        parts.append('cycle number %d' % cycle)

    return _described(self, ', '.join(parts))


# ----- observer events ---------------------------------------------------------------

def _where_the_tracker_is(event):
    """`at t=0.31 in 30 digits`, for any event that carries a tracker."""
    tracker = _ask(lambda: event.tracker())
    if tracker is None:
        return None

    point = _ask(tracker.current_point)
    if point is None or len(point) == 0:
        return None

    where = 'at t=%s' % _time(_ask(tracker.current_time))
    digits = _ask(tracker.current_precision)
    return '%s in %d digits' % (where, digits) if digits else where


def _event_repr(self):
    return _described(self, _where_the_tracker_is(self))


def _precision_event_repr(self):
    parts = []
    previous, following = _ask(lambda: self.previous()), _ask(lambda: self.next())
    if previous is not None and following is not None:
        parts.append('%d -> %d digits' % (previous, following))

    tracker = _ask(lambda: self.tracker())
    if tracker is not None:
        point = _ask(tracker.current_point)
        if point is not None and len(point):
            parts.append('at t=%s' % _time(_ask(tracker.current_time)))

    return _described(self, ', '.join(parts))


def _stepsize_event_repr(self):
    parts = []
    tracker = _ask(lambda: self.tracker())
    if tracker is not None:
        stepsize = _ask(tracker.current_stepsize)
        if stepsize is not None:
            parts.append('step now %s' % _number(stepsize))
    where = _where_the_tracker_is(self)
    if where:
        parts.append(where)
    return _described(self, ', '.join(parts))


def _path_event_repr(self):
    parts = []
    index = _ask(lambda: self.path_index())
    if index is not None:
        parts.append('path %d' % index)
    where = _where_the_tracker_is(self)
    if where:
        parts.append(where)
    return _described(self, ', '.join(parts))


# ----- collectors --------------------------------------------------------------------

def _path_data_collector_repr(self):
    times = _ask(self.times)
    steps = 0 if times is None else len(times)
    return _described(self, _count(steps, 'step') if steps else 'nothing recorded')


def _solution_path_collector_repr(self):
    series = getattr(self, 'series', None)
    if not series:
        return _described(self, 'nothing recorded')

    steps = 0
    for path in series:
        recorded = _ask(path.times)
        steps += 0 if recorded is None else len(recorded)
    return _described(self, '%s, %s' % (_count(len(series), 'path'), _count(steps, 'step')))


def _sample_collector_repr(self):
    runs, samples = _ask(self.num_runs), _ask(self.num_samples)
    if not runs and not samples:
        return _described(self, 'nothing recorded')
    return _described(self, '%s, %s' % (_count(runs, 'run', 'no runs'),
                                        _count(samples, 'sample', 'no samples')))


# ----- linear algebra ----------------------------------------------------------------

def _decomposition_repr(self):
    rows, cols = _ask(self.rows), _ask(self.cols)
    if rows is None or cols is None:
        return _described(self, None)
    return _described(self, '%dx%d' % (rows, cols))


# ----- the rest ----------------------------------------------------------------------

def _midpath_report_repr(self):
    passed = _ask(lambda: self.passed())
    crossings = _ask(lambda: self.num_crossings_detected())
    attempts = _ask(lambda: self.num_resolve_attempts())

    parts = ['passed' if passed else 'did not pass']
    if crossings is not None:
        parts.append(_count(crossings, 'crossing'))
    if attempts:
        parts.append(_count(attempts, 're-track'))
    return _described(self, ', '.join(parts))


def _user_start_system_repr(self):
    return _described(self, _count(_ask(lambda: self.num_start_points()), 'start point'))


def _bare_repr(self):
    """For the types that expose nothing to say anything about."""
    return _described(self, None)


# ----- containers of nodes -----------------------------------------------------------

def _sequence_repr(self):
    """A list of the elements' own reprs.

    The bound containers stream their elements in C++, and an element that is a pointer
    to a node streams as its address -- so a variable group read as `[0x55f1a4, 0x55f1b0]`
    while its str() read `[x,y]`.
    """
    try:
        inside = ', '.join(repr(element) for element in self)
    except Exception:
        return '<%s>' % type(self).__name__
    return '%s([%s])' % (type(self).__name__, inside)


def _plain_sequence_repr(self):
    try:
        return '[%s]' % ', '.join(repr(element) for element in self)
    except Exception:
        return '<%s>' % type(self).__name__


# ----- installation ------------------------------------------------------------------

def install():
    """Attach every repr.  Idempotent, and safe to call before or after anything else."""
    import bertini._pybertini as native

    _install_trackers(native)
    _install_endgames()
    _install_events(native)
    _install_observers(native)
    _install_collectors(native)
    _install_linalg(native)
    _install_containers(native)
    _install_remainder(native)


def _install_trackers(native):
    _each(_classes(native.tracking, ('AMPTracker', 'DoublePrecisionTracker',
                                     'MultiplePrecisionTracker')), _tracker_repr)


def _install_endgames():
    import bertini.endgame as endgame

    # the python wrappers forward every accessor to the C++ endgame they hold, so one
    # function on the shared base serves all six
    base = getattr(endgame, '_EndgameBase', None)
    _attach(base, _endgame_repr)
    _each(_classes(endgame, ('AMPCauchyEndgame', 'AMPPowerSeriesEndgame',
                             'FixedDoubleCauchyEndgame', 'FixedDoublePowerSeriesEndgame',
                             'FixedMultipleCauchyEndgame', 'FixedMultiplePowerSeriesEndgame')),
          _endgame_repr)


_TRACKING_FLAVOURS = ('amp', 'double', 'multiple')


def _install_events(native):
    for flavour in _TRACKING_FLAVOURS:
        module = getattr(native.tracking.observers, flavour, None)
        if module is None:
            continue

        _each(_classes(module, ('PrecisionChanged', 'PrecisionIncreased',
                                'PrecisionDecreased')), _precision_event_repr)
        _each(_classes(module, ('StepsizeIncreased', 'StepsizeDecreased')),
              _stepsize_event_repr)
        _each(_classes(module, ('TrackingEvent', 'TrackingStarted', 'TrackingEnded',
                                'SuccessfulStep', 'FailedStep', 'InfinitePathTruncation')),
              _event_repr)

    algorithm = getattr(native.nag_algorithms, 'observers', None)
    if algorithm is not None:
        _each(_classes(algorithm, ('PathStarted', 'PathComplete')), _path_event_repr)
        _each(_classes(algorithm, ('AlgorithmEvent', 'AlgorithmStarted',
                                   'AlgorithmComplete')), _bare_repr)

    # the event a custom observer is handed when nothing more specific fits.  It is not
    # exported under a name of its own, so it is only reachable here.
    _attach(getattr(native.detail, 'AnyEvent', None), _bare_repr)


def _install_observers(native):
    for flavour in _TRACKING_FLAVOURS:
        module = getattr(native.tracking.observers, flavour, None)
        if module is not None:
            _each(_classes(module, ('CustomObserver', 'GoryDetailLogger',
                                    'FirstPrecisionRecorder')), _bare_repr)

    for flavour in ('amp_cauchy', 'amp_pseg', 'double_cauchy', 'double_pseg',
                    'multiple_cauchy', 'multiple_pseg'):
        module = getattr(native.endgame.observers, flavour, None)
        if module is not None:
            _each(_classes(module, ('CustomObserver', 'GoryDetailLogger')), _bare_repr)

    algorithm = getattr(native.nag_algorithms, 'observers', None)
    if algorithm is not None:
        _each(_classes(algorithm, ('CustomObserver',)), _bare_repr)

    import bertini.tracking as tracking
    _each(_classes(tracking, ('CallbackObserver', 'PathCollectionObserver')), _bare_repr)
    _attach(getattr(native.tracking, 'Observable', None), _bare_repr)


def _install_collectors(native):
    import bertini.tracking as tracking
    import bertini.nag_algorithm as nag

    _attach(getattr(tracking, 'PathDataCollector', None), _path_data_collector_repr)
    for flavour in _TRACKING_FLAVOURS:
        module = getattr(native.tracking.observers, flavour, None)
        if module is not None:
            _each(_classes(module, ('PathDataCollector',)), _path_data_collector_repr)

    _attach(getattr(nag, 'SolutionPathCollector', None), _solution_path_collector_repr)

    for flavour in ('amp_cauchy', 'amp_pseg', 'double_cauchy', 'double_pseg',
                    'multiple_cauchy', 'multiple_pseg'):
        module = getattr(native.endgame.observers, flavour, None)
        if module is not None:
            _each(_classes(module, ('SampleSequenceCollector',)), _sample_collector_repr)


def _install_linalg(native):
    _each(_classes(native.linalg, ('ColPivHouseholderQR', 'ColPivHouseholderQRReal',
                                   'HouseholderQR', 'HouseholderQRReal',
                                   'JacobiSVD', 'JacobiSVDReal',
                                   'PartialPivLU', 'PartialPivLUReal')),
          _decomposition_repr)


def _install_containers(native):
    # these have a repr already; it is the one being fixed, since it streams the elements
    # in C++ and an element that is a pointer to a node streams as its address
    _attach(getattr(native.container, 'VariableGroup', None), _sequence_repr, replacing=True)
    _each(_classes(native.container, ('ListOfVariableGroup', 'ListOfComplex',
                                      'ListOfNode', 'ListOfRational')),
          _plain_sequence_repr, replacing=True)


def _install_remainder(native):
    _attach(getattr(native.nag_algorithms, 'MidpathCheckReport', None), _midpath_report_repr)
    _attach(getattr(native.nag_algorithms, 'UserStartSystem', None), _user_start_system_repr)
    _each(_classes(native.nag_algorithms, ('StartSystemFactory', 'AnyZeroDim')), _bare_repr)

    # the decomposition solvers are scaffolding whose Solve() throws; there is nothing to
    # report about one yet, but an address is still the wrong thing to print
    _each([getattr(native.nag_algorithms, name) for name in dir(native.nag_algorithms)
           if name.startswith('NID')
           or name.startswith('NumericalIrreducibleDecomposition')], _bare_repr)
