=============================================================================================
Tools for MultiXrank
=============================================================================================

MultiXrank is a Python package dedicated to the exploration of heterogeneous multilayer networks (aka universal multilayer networks) with random walk with restart. MultiXrank can be used in a wide range of applications. The Multixrank package is accessible there: https://github.com/anthbapt/multixrank
If you use MultiXrank, **please cite the following article**:

**Baptista, A., González, A., Baudot, A.**.
Baptista, A., Gonzalez, A. & Baudot, A. Universal multilayer network exploration by random walk with restart. Commun Phys 5, 170 (2022)

https://doi.org/10.1038/s42005-022-00937-9

You can find here the scripts used in the paper and adaptable for your own applications of MultiXrank:

.. code-block:: bash
    :caption:
        1: Scripts for the evaluation (**evaluation**)
        .. hlist::
            * Leave-One-Out cross validation
            * Link prediction

        2: Scripts for the parameter exploration (**exploration**)
        
        3: Further scripts will be added soon (07/11/2022)

The scripts use the multiprocessing python library. You can specify the number of threads dedicated to the computation.
