from .utils import *
from .save_uniem_annotator import *

from .BuildAnyShape import *
from .CreateVolumeFromReconstruct import *
from .RotateVolume import *
from .GetLabeledConcs import *
from .ConnectLabeledConcs import *
from .ConnectAnalysis import *
from .ConnectTotalConcs import *
from .CreateLabelVolumeFromUniEM import *
from .CreateSurface import *
from .LoadLabelVolume import *
from .RepeatRun import *

"""
__all__ = ['buildAnyShape','calcSurfaceArea','createVolumeFromReconstruct','rotateVolume']

from . import utils
from inspect import getmembers, isfunction
function_list_utils = [o[0] for o in getmembers(utils) if isfunction(o[1])]
__all__.extend(function_list_utils)
"""


import os, glob, importlib
from inspect import getmembers, isfunction, isclass

base_dir     = os.path.dirname(__file__)
python_files = glob.glob(os.path.join(base_dir, '[a-zA-Z0-9]*.py'))
targ_names   = [os.path.split(os.path.splitext(file)[0])[1] for file in python_files ]

__all__ = []
for targ_name in targ_names:
	m = importlib.import_module('pyLD.'+targ_name)
	classes_funcs = [o[0] for o in getmembers(m) if isfunction(o[1]) or isclass(o[1])]
	classes_funcs = [s for s in classes_funcs if s not in ['flatten', 'rescale']]
	__all__.extend(classes_funcs)
	# print(', '.join(classes_funcs) )

print('Importing pyLD files:')
print(', '.join(targ_names) )

# https://stackoverflow.com/questions/46263274/how-can-i-prevent-sphinx-from-displaying-the-full-path-to-my-class
# https://qiita.com/suzuki-kei/items/8fea67655abf216a5013


