🎲 Continuing a computation with the same random draws
*******************************************************

.. testsetup:: *

   import pickle
   import bertini

A long computation runs in stages, and every stage draws random numbers: the gamma of a homotopy, a patch, the coefficients of a start system, a random projection. If the computation writes itself out between stages, so that it can stop and continue later, it must also write out where its random draws had got to. Otherwise the stages after a restart draw different numbers than an uninterrupted run would have, and the continued run is a different computation.

A seed is not enough for this: it says where the draws *began*, not where they have reached. The **random state** says where they have reached.

Getting and setting the random state
====================================

:func:`bertini.random.get_random_state` returns the calling thread's random state, and :func:`bertini.random.set_random_state` sets it. The draws after setting a state are the draws that followed getting it:

.. testcode::

   bertini.random.set_random_seed(1042)
   bertini.random.random_vector(3)          # the first stage's draws

   state = bertini.random.get_random_state()
   after_the_first_stage = list(bertini.random.random_vector(3))

   bertini.random.set_random_seed(7)        # something else draws in between
   bertini.random.random_vector(5)

   bertini.random.set_random_state(state)
   assert list(bertini.random.random_vector(3)) == after_the_first_stage

Solves draw from the same state, so a solve after setting a state builds the same homotopy, and finds the same solutions, as the solve that followed getting it.

What the state holds
====================

A :class:`~bertini.random.RandomState` holds everything the thread's next draws depend on: the global seed, whether one has been chosen, the ordinal of the next solve seed (:func:`~bertini.random.derive_solve_seed`), and the thread's draw stream. Two states are equal when all of these are:

.. testcode::

   bertini.random.set_random_seed(5)
   bertini.random.derive_solve_seed()
   state = bertini.random.get_random_state()
   print(state.global_seed, state.seed_chosen, state.solve_ordinal)
   assert bertini.random.get_random_state() == state

.. testoutput::

   5 True 1

The fields can be set, so a state can be built field by field. Each path of a solve, and each worker of a parallel solve, draws from its own stream, reseeded from the global seed and its index when it starts, so those streams are not part of the state.

Writing the state out
=====================

A state is written as versioned plain text, and read back from it, so it can sit beside the rest of a computation's output:

.. testcode::

   text = state.to_text()
   print(text.splitlines()[0])
   assert bertini.random.RandomState.from_text(text) == state

.. testoutput::

   b2randstate/1

A state also pickles (as its text), so it can go wherever a pickle goes:

.. testcode::

   assert pickle.loads(pickle.dumps(state)) == state

Text of another version, or text that is not a state, is refused:

.. testcode::

   try:
       bertini.random.RandomState.from_text("not a random state")
   except ValueError as e:
       print("refused")

.. testoutput::

   refused

In C++, the same state is ``bertini::RandomState``, with ``GetRandomState()``, ``SetRandomState(state)``, ``ToText()`` and ``RandomState::FromText(text)``; it also goes through Boost.Serialization archives, and so through MPI.
