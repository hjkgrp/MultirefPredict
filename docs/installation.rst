Installation
============

Currently MultirefPredict can only be installed from source code

1. Prerequisite: have `Anaconda or miniconda <https://www.anaconda.com/distribution/>`_ installed on your system

2. Clone MultirefPredict source from gibhub

   .. code-block:: shell

      git clone https://github.com/hjkgrp/MultirefPredict.git

3. Go to the folder root folder for MultirefPredict, create the conda environment from yaml file

   .. code-block:: shell
   
      cd dev_tools/conda_envs
      conda env create -f test_env.yaml

   .. IMPORTANT::
      It is highly recommended to install the packages needed by creating environment from this yaml file. This help
      make sure that the version of package installed (e.g. psi4) is the same as the one tested by the developer.
      Installing the required packages from other distributions may result in unexpected behavior

4. Now you have created an environment called **test** for with the prerequisite packages for running MultirefPredict. 
   You can rename the environment as you like by cloning this environment to a new name and remove the origional one.
   For example, to rename it to **multiref**:

   .. code-block:: shell

      conda create --name multiref --clone test
      conda remove --name test --all

5. Activate the conda environment you just created. Go back to the root directory of MultirefPredict (where the setup.py
   file locates). Local install with pip

   .. code-block:: shell
      
      cd Your_root_path_to_MultirefPrect
      pip install -e .

