👋 Welcome to Bertini 2
====================================

Bertini 2 is software for numerically solving systems of polynomials.

🗺 Mathematical overview
----------------------------------

The main algorithm for numerical algebraic geometry implemented in Bertini is homotopy continuation.  A homotopy is formed, and the solutions to the start system are continued into the solutions for the target system.


.. figure:: images_common/homotopycontinuation_generic.png
   :scale: 100 %
   :height:  285 px
   :width: 400 px
   :alt: Homotopy continuation

   Predictor-corrector methods with optional adaptive precision track paths from 1 to 0, solving :math:`f`.

The definitive resource for Bertini 1 is the book :cite:`bertinibook`.  While the way we interact with Bertini changes from version 1 to version 2, particularly when using the Python library, the algorithms remain fundamentally the same.  So do most of the ways to change settings for the path trackers, etc.  We believe that embracing the flexibility of Python allows for much greater flexibility.  It also will relieve the user from the burden of input and output file writing and parsing.  Instead, computed results are returned directly to the user.

Consider checking out the :ref:`🔦 Tutorials <tutorials>`.

📥 Installation
----------------

Easiest is to install from `pip`:

::

    pip install bertini2

Note that the import in Python does NOT use the 2, as silviana didn't want to have that appear everywhere in the code.  So it's just

::

    import bertini

and then do whatever with the library.

.. note::

   You can find more complete instructions, including instructions around MPI parallelism, in `the Wiki <https://github.com/bertiniteam/b2/wiki/Installation>`_.



⛲️ Source code
---------------------------

The Bertini 2 source code is available at `its GitHub repo <https://github.com/bertiniteam/b2>`_.

The core is written in template-heavy C++, and is exposed to Python through Boost.Python.

⚖️ Licenses
------------------

Bertini2 and its direct components are available under GPL3, with additional clauses in section 7 to protect the Bertini name.  Bertini2 also uses open source softwares, with their own licenses, which may be found in the Bertini2 repository, in the licenses folder.

🤖 Notes on use of AI
--------------------------------------

Silviana Amethyst has been using AI Coding agents since May 2026 to assist in the writing of Bertini 2.  I make no effort to obfuscate my use of these tools.  I am using them, and it's been amazing.

Here are some ways I am doing so:

* Bug fixing.  I have fixed dozens and dozens of bugs, from subtle to obvious, using these tools.
* Performance analysis.  Claude wrote a LU factorization replacement in about 30 minutes, replacing the one from Eigen, and it improved the performance of the mpfr path tracking significantly.
* Help writing documentation.  I used Claude to write first versions of some of the tutorials, and I continue to use it to revise those tutorials as interface changes.  I have hand-edited the tutorials to be closer to my own style, and I have scientific faith in them.
* Continuous integration.  I was trapped in a situation, where almost no one could use this codebase because they would never succeed in installing it.  After Hong Kee at CBG helped me get over the hump so I had basics of CI running, and observing the ways he used OpenCode to help, I took a dive.  I was able to solve a decade-old barrier.  Now I have full CI running on every target OS i think is reasonable.  I have automatic documentation generation.  I have checks that help improve the quality of the code before it lands in a release.

My use of AI coding tools has hinged fully on my expertise in numerical algebraic geometry.  I've been working on software in this real since at least 2009.  Engaging with AI on this project has helped propel me towards the success I've always dreamed of.

There are some parts of this library that I extensively used Claude Fable 5 to assist with design: hashconsing of the symbolic math part of the library, and with the records system.  I have had dreams for many years of these two things.  AI made it possible for me.  I mostly work on this set of tools in my evenings and weekends.   I would never have had the time to accomplish these goals without these tools.  I hope you get some benefit from them, too.

I assume responsibility for the mathematical and computational correctness of this library.  I assume the credit for the work in my commits where Claude or another agentic system is listed.  I keep detailed records of my work, particularly my prompts.  I remain dedicated to mathematical and scientific excellence, both with and without the use of AI coding tools.
