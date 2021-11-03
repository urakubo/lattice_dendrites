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

  # System paths
  export LD_LIBRARY_PATH=/lib64:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH

  # Requirements in $LM_INSTALLED/hdf5_1.12.0_gcc8.4.0/lib:
  # libhdf5.so.200, libhdf5_hl.so.200
  # Requirements in /usr/local/cuda/lib64:
  # libcudart.so.11.0
  # Requirements in /lib64:
  # ld-linux-x86-64.so.2
  # Requirements in /usr/lib/x86_64-linux-gnu:
  # libstdc++.so.6, libpython3.6m.so.1.0
  # Requirements in /lib/x86_64-linux-gnu:
  # libm.so.6, libgcc_s.so.1, libpthread.so.0, libc.so.6, libz.so.1, libdl.so.2, librt.so.1, libexpat.so.1, libutil.so.1

```

This is a case of bash. You can also write the above setting to .bashrc .

