Running CAMPAREE
================

.. _running-cmdoptions:


Command Line Options
--------------------

.. code-block:: none

    usage: run_camparee.py [-h] -c CONFIG [-r RUN_ID] [-d] [-m {lsf,serial,sge}]
                           [-s SEED]

    CAMPAREE - RNA molecule simulator (v0.2.0)

    optional arguments:
      -h, --help            show this help message and exit

    required named arguments:
      -c CONFIG, --config CONFIG
                            Full path to configuration file.

    optional named arguments - these override configuration file arguments.:
      -r RUN_ID, --run_id RUN_ID
                            Alphanumberic used to specify run id (letters, numbers
                            and underscores only).
      -d, --debug           Indicates whether additional diagnostics are printed.
      -m {lsf,serial,sge}, --scheduler_mode {lsf,serial,sge}
                            Indicates whether to dispatch jobs serially, or using
                            a job scheduler
      -s SEED, --seed SEED  Optional integer value used as a seed for random
                            number generation.
