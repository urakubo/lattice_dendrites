import glob
import os
from os.path import join as pjoin
import subprocess
import sys

def get_cuda_version(cuda_home):
    """Locate the CUDA version
    """
    if not isinstance(cuda_home, str):
        return None
    version_file = os.path.join(cuda_home, "version.txt")
    try:
        if os.path.isfile(version_file):
            with open(version_file) as f:
                version_str = f.readline().replace('\n', '').replace('\r', '')
                return version_str.split(" ")[2][:4]
        else:
            version_str = subprocess.check_output([os.path.join(cuda_home,"bin","nvcc"),"--version"])
            version_str=str(version_str).replace('\n', '').replace('\r', '')
            idx=version_str.find("release")
            return version_str[idx+len("release "):idx+len("release ")+4]
    except:
        raise RuntimeError("Cannot read cuda version file") 
def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'include' and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDA_HOME or CUDA_PATH env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """
    # Guess #1
    cuda_home = os.environ.get('CUDA_HOME') or os.environ.get('CUDA_PATH')
    if cuda_home is None:
        # Guess #2
        try:
            nvcc = subprocess.check_output(
                ['which', 'nvcc']).decode().rstrip('\r\n')
            cuda_home = os.path.dirname(os.path.dirname(nvcc))
        except subprocess.CalledProcessError:
            # Guess #3
            cuda_home = '/usr/local/cuda'
            if not os.path.exists(cuda_home):
                cuda_home = None
    version = get_cuda_version(cuda_home)

    return version

if __name__ == '__main__':

    print('CUDA_VERSION: ', locate_cuda())
