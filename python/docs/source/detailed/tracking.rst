🛤 bertini.tracking
===========================


Notes
--------

Trackers in Bertini2 are stateful objects, that refer to a system they are tracking, hold their specific settings, and have a notion of current time and space value.

Here are some particular classes and functions to pay attention to:

* :class:`bertini.tracking.AMPTracker`
* :class:`bertini.tracking.DoublePrecisionTracker`
* :class:`bertini.tracking.MultiplePrecisionTracker`

Here are the implemented ODE predictors you can choose from:

* :class:`bertini.tracking.Predictor`

Calls to :meth:`track_path` return a :class:`bertini.tracking.SuccessCode`.

And, trackers are implemented using observer pattern.  They live in the ``bertini.tracking.observers`` namespace, with provisions for each tracker type available under a submodule thereof: ``amp``, ``multiple``, and ``double``.  They are also conveniently available using the ``tr.observers``, where ``tr`` is a tracker you already made.  See :mod:`bertini.tracking.observers.amp`


Auto-generated docs
--------------------

.. automodule:: bertini.tracking


🛤 bertini.tracking.observers
===================================

📝 All of these are available for all trackers, though you should use the ones for your tracker type.  Look in ``bertini.tracking.AMPTracker.observers``, etc.

.. automodule:: bertini.tracking.observers

#. ``bertini.tracking.observers.amp``
#. ``bertini.tracking.observers.double``
#. ``bertini.tracking.observers.multiple``

📝 Symmetrically, there are the same observers in all three.

.. automodule:: bertini.tracking.observers.amp


 Know that you are loved and appreciated, dear reader.  💟
