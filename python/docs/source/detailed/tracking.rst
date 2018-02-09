🛤 pybertini.tracking
====================

.. include:: common_doc_nav.incl

Notes
--------

Trackers in Bertini2 are stateful objects, that refer to a system they are tracking, hold their specific settings, and have a notion of current time and space value.

Here are some particular classes and functions to pay attention to:

* :class:`pybertini.tracking.AMPTracker`
* :class:`pybertini.tracking.DoublePrecisionTracker`
* :class:`pybertini.tracking.MultiplePrecisionTracker`

Here are the implemented ODE predictors you can choose from:

* :class:`pybertini.tracking.Predictor`

Calls to :meth:`track_path` return a :class:`pybertini.tracking.SuccessCode`.

And, trackers are implemented using observer pattern.  They live in the ``pybertini.tracking.observers`` namespace, with provisions for each tracker type available under a submodule thereof: ``amp``, ``multiple``, and ``double``.  They are also conveniently available using the ``tr.observers``, where ``tr`` is a tracker you already made.  See :mod:`pybertini.tracking.observers.amp`


Auto-generated docs
--------------------

.. automodule:: pybertini.tracking


🛤 pybertini.tracking.config
==========================

.. automodule:: pybertini.tracking.config


🛤 pybertini.tracking.observers
================================

📝 All of these are available for all trackers, though you should use the ones for your tracker type.  Look in ``pybertini.tracking.AMPTracker.observers``, etc.

.. automodule:: pybertini.tracking.observers

#. ``pybertini.tracking.observers.amp``
#. ``pybertini.tracking.observers.double``
#. ``pybertini.tracking.observers.multiple``

📝 Symmetrically, there are the same observers in all three.

.. automodule:: pybertini.tracking.observers.amp


 Know that you are loved and appreciated, dear reader.  💟
