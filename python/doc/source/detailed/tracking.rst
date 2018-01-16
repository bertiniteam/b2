pybertini.tracking
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

Auto-generated docs
--------------------

.. automodule:: pybertini.tracking


pybertini.tracking.config
==========================

.. automodule:: pybertini.tracking.config
