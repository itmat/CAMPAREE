Running CAMPAREE
================

Installation
------------

Installation Prerequisites:

- git
- python version 3.6
- Java 1.8

The following installation instructions are for Python3 installed in a Linux
environment. Windows/Mac environments, or environments using an Anaconda/Miniconda
installation of Python may require slight modifications to the commands below.

1. Download CAMPAREE and BEERS_UTILS from their respective git repos::

    git clone https://github.com/itmat/CAMPAREE.git
    git clone https://github.com/itmat/BEERS_UTILS.git

2. Create a Python virtual environment to install required Python libraries to::

    cd CAMPAREE
    python3 -m venv ./venv_camparee

3. Activate the Python environment::

    source ./venv_camparee/bin/activate

4. Install required libraries::

    pip install -r requirements.txt
    pip install -r ../BEERS_UTILS/requirements.txt

5. Install CAMPAREE and BEERS_UTILS packages in the Python environment::

    pip install -e .
    pip install -e ../BEERS_UTILS

6. Prepare/install resource files for the organism of choice (see full instructions in *Resource Files* section)


Command Line Options
--------------------

.. code-block:: none

    usage: run_camparee.py [-h] -c CONFIG [-r RUN_ID] [-d] [-m {lsf,serial,sge}]
                           [-s SEED]

    CAMPAREE - RNA molecule simulator (v0.1.0)

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

Running CAMPAREE on AWS
-----------------------

While CAMPAREE can run on a single machine, it can require 50-60 GB of RAM to
simulated samples from a mammalian-size genome, and a cluster computing environment
can substantially decrease run times. For users who do not have access to machines
with enough memory or a supported cluster environment, Amazon Web Services (AWS)
provide a viable alternative for running CAMPAREE. These instructions assume users
have access to an AWS account with associated key pairs (i.e. \*.pem file), and
have a passing familiarity with checking out AWS instances through the web
interface (the AWS Management Console).

Note: AWS is an excellent resource for providing pay-as-you-go access to high-end
machines. However, as a cloud-based service, it is important users take approrpiate
security measure to secure their AWS instances. A full discussion of the many
security feature available through AWS, such as restricting access to a range of
IP addresses, is the beyond the scope of this document. These instructions are
meant as a starting point, and are written to be as broadly applicable as possible.
For example, the instructions below create instances with publically accessible
IP addresses and AWS-managed security keys. This is likely sufficient security for
users attempting to run the simulator with published data.

Here are instructions for two possible scenarios for running CAMPAREE on AWS.

Single Instance in Serial Mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By checking out a single, high-memory machine, users can run CAMPAREE on AWS in
serial mode. Note, the instructions listed here for checking out an AWS instance
rely mostly on the default settings. Users with managed AWS accounts may need
to alter their accounts to allow them to checkout out instances, and/or adapt
these instructions to work within their accounts restrictions.

As RAM tends to be the bottle neck for most CAMPAREE processes, rather than the
number of CPU cores, we recommend users select a `Memory Optimized <https://aws.amazon.com/ec2/instance-types/>`_
machine, like the r4.2xlarge (8 vCPUs; 61 GiB RAM).

1. Configure and checkout the AWS instance

    AMI:
        ami-07ebfd5b3428b6f4d (64-bit x86) - Ubuntu Server 18.04 LTS (HVM), SSD Volume Type

    Configure Instance Details (defaults except for the following):
        - Network: Depends on user's account.
        - Subnet: Depends on user's account. Must be in same region as the desired key pair.
        - Auto-assign Public IP: Enable

    Add Storage:
        Default volume:
            Size: *Still working out recommendations based on input size*

    Select Key pair:
        Choose an existing key pair (will only display options from the same regions as the subnet selected above).

2. Confirm instance details and wait for instance to enter running state.

3. Login to instance with ssh (requires the \*.pem file associated with the key selected above).

4. Install CAMPAREE pre-requisites::

    sudo apt-get update && sudo apt-get -y upgrade
    sudo apt-get -y install openjdk-8-jre python3-venv python3-pip

5. Follow CAMPAREE installation instructions listed above. Note, you can test the install using the 'Baby Genome' as described in the quick start guide.

6. Prepare CAMPAREE config file, making sure to set ``scheduler_mode:`` to 'serial'.

7. CAMPAREE is now ready to run in **serial** mode.

AWS ParallelCluster in SGE Mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `AWS ParallelCluster <https://aws.amazon.com/hpc/parallelcluster/>`_ program
provides a command line utility to create and manage a cluster environment entirely
on AWS. ParalleleCluster will automatically add and remove compute nodes as they
are needed, and supports several job managers (including SGE). Currently, the
account used to launch the ParallelCluster must have full admin privileges (the
default for most unmanaged AWS accounts).

Again, we recommend using `Memory Optimized <https://aws.amazon.com/ec2/instance-types/>`_
machines, like the r4.2xlarge, for the compute nodes on the ParallelCluster. For
the master node, a machine with lower specs, like the m5.large (2 vCPUs; 8 GiB
RAM), should be adequate.

1. Install AWS ParallelCluster and its prerequisites according to `these instructions <https://docs.aws.amazon.com/parallelcluster/latest/ug/install.html>`_.

2. Configure ParallelCluster by running the ``pcluster configure`` command and entering the following options

    - AWS Region ID: <enter region matching desired AWS key pair>
    - EC2 Key Pair Name: <select desired AWS key pair>
    - Scheduler: 'sge'
    - Operating System: 'ubuntu1804' (required for pre-install script used below)
    - Minimum cluster size (instances): '1'
    - Maximum cluster size (instances): '10' (By default, extra instances are only added when a job waits in the queue for 10 minutes)
    - Master instance type: m5.large
    - Compute instance type: r4.2xlarge
    - Automate VPC creation? y
    - Network Configuration: Master in a public subnet and compute fleet in a private subnet
    - Automate Subnet creation? y

3. Add custom startup script to the ParallelCluster config.

    ParallelCluster supports post-install scripts that run on each of the nodes
    after they've been added to the cluster. This script will install all of
    CAMPAREE's prerequisites. Start by opening the ParallelCluster config file
    (generally located at ``~/.parallelcluster/config``) with a text editor.
    Append the following lines to the end of the ``[cluster default]`` section
    of the config file::

            post_install = s3://itmat.data-simulators/parallelcluster_camparee_prereqs_postinstall_Ubuntu1804.sh

4. Launch the ParallelCluster::

    pcluster create camparee-cluster

5. Once the cluster is full initialized, connect to the master node::

    pcluster ssh camparee-cluster -i /path/to/AWS_key_file.pem

6. Install CAMPAREE on the cluster using the instructions listed above. Note, the prerequisites were already handled by the post-install script. You can test the install using the 'Baby Genome' as described in the quick start guide.

7. Prepare CAMPAREE config file, making sure to set ``scheduler_mode:`` to 'sge'.

8. CAMPAREE is now ready to run in **sge** mode.

When you have finish running CAMPAREE and have transferred all data off of the
cluster, you can shut down and delete the cluster with the following command::

    pcluster delete camparee-cluster

Note, all data on the ParallelCluster will be lost after this command completes.
