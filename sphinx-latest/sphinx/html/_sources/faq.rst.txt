Troubleshooting and FAQ
========================

.. contents::
   :local:
   :depth: 1


Troubleshooting
----------------

It is taking Conda forever to solve the environment during the GPlately installation.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If Conda gets **stuck while solving the environment** during the installation of GPlately, you can try to use micromamba_ instead.

Unable to install or use GPlately on my computer.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For various reasons, GPlately may not install or function correctly on some computers. 
To address this, we've prepared a Docker image with a fully working GPlately installation. 
See `this web page <installation.html#use-docker>`__.

Failed to import PyGMT.
~~~~~~~~~~~~~~~~~~~~~~~

PyGMT requires Python>=3.11. Upgrade your Python if you want PyGMT.

Parallel code running slower than serial code.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is because the overhead of multiprocessing is too high, outweighing the benefits of parallelization.
Optimizing your workflow for parallel computing can significantly improve performance. For example, 
pickling certain Python objects can be costly. In some cases, it's faster to recreate these objects than to serialize them.

Frequently Asked Questions
--------------------------

Where to get help?
~~~~~~~~~~~~~~~~~~

You may post questions on the `GPlates online forum <https://discourse.gplates.org/>`__.

Where to report bugs?
~~~~~~~~~~~~~~~~~~~~~

You may `create issues <https://github.com/GPlates/gplately/issues>`__ in the GitHub GPlately repository.

How to contact the development team?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You may `contact EarthByte research group <https://www.earthbyte.org/contact-us-3/>`__.

Is it safe to use GPlately?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Best effort has been made to make the software safe. However, the GPlately software is provided "as is", without warranty.
If you have legal or security concerns, you may `contact EarthByte research group <https://www.earthbyte.org/contact-us-3/>`__.


.. _micromamba: https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html