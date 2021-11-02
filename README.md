# lattice_dendrites
Lattice Dendrites (LD) is an extension of lattice microbes (LM) that can incorporate morphologically realistic shapes of cells, in particular, from 3-dimentional (3D) reconstruction of electron microscopic (EM) images. LD and LM work on Python2. I will rewrite them to work on Python3. For detailed information, please see the [manual](https://urakubo.github.io/lattice_dendrites/).


## Installation

We parepared LM for py3.6/CUDA11.0/ubuntu18.04, you can retrieve it by

  > tar -xvzf lm_py3.6_cuda11.0_ubuntu18.04.tar.gz

Then, create paths as follows:

```bash
  # LM paths
  export PATH=$LM_INSTALLED/bin:$PATH
  export PYTHONPATH=$LM_INSTALLED/lib/lm:$PYTHONPATH
  export PYTHONPATH=$LM_INSTALLED/lib/python:$PYTHONPATH
  export LD_LIBRARY_PATH=$LM_INSTALLED/hdf5_1.12.0_gcc8.4.0/lib:$LD_LIBRARY_PATH

  # CUDA paths
  export PATH=/usr/local/cuda/bin:$PATH
  export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

  # System paths (unnecessay in many cases)
  export LD_LIBRARY_PATH=/lib64:/usr/x86_64-linux-gnu:/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH
```

This is a case of bash. You can also write the above setting to .bashrc .

