from .Utils import *
from .PlotCompartmentModelBackend import PlotCompartmentModelBackend, Interactor
from .SimulateMembranePotential import SimulateMembranePotential
from .CreateGraph import CreateGraph
from .CloseMesh import CloseMesh
from .DendriticCompartments  import DendriticCompartments
from .SpineCompartments      import SpineCompartments
from PlotCompartmentModelBackend import PlotCompartmentModelBackend, Interactor


#'''
import os, glob, importlib
from inspect import getmembers, isfunction, isclass

base_dir     = os.path.dirname(__file__)
python_files = glob.glob(os.path.join(base_dir, '[a-zA-Z0-9]*.py'))
targ_names   = [os.path.split(os.path.splitext(file)[0])[1] for file in python_files ]
targ_names   = [name for name in targ_names if name not in  ['CloseMesh','PlotCompartmentModelBackend'] ]

print('targ_names ', targ_names)

__all__ = []
for targ_name in targ_names:
	m = importlib.import_module('pyLD.'+targ_name)
	classes_funcs = [o[0] for o in getmembers(m) if isfunction(o[1]) or isclass(o[1])]
	classes_funcs = [s for s in classes_funcs if s not in ['flatten', 'rescale']]
	__all__.extend(classes_funcs)
	# print(', '.join(classes_funcs) )

print('Importing pyLD files:')
print(', '.join(targ_names) )
#'''

# https://stackoverflow.com/questions/46263274/how-can-i-prevent-sphinx-from-displaying-the-full-path-to-my-class
# https://qiita.com/suzuki-kei/items/8fea67655abf216a5013


