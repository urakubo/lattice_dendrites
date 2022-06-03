from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install

import os, sys
import tarfile
from pyLD.check_cuda import locate_cuda


with open("README.md") as f:
	long_description = f.read()

class PostInstallCommand(install):
	"""Post-installation for installation mode."""
	def run(self):
		install.run(self)
		# PUT YOUR POST-INSTALL SCRIPT HERE or CALL A FUNCTION
		# print(self.install_lib)


class PostDevelopCommand(develop):
	"""Post-installation for development mode."""
	def run(self):
		develop.run(self)
		# PUT YOUR POST-INSTALL SCRIPT HERE or CALL A FUNCTION

s = setup(
	name="lattice_dendrites",
	version="0.1.0",
	author="hidetoshi-urakubo",
	author_email="hurakubo@gmail.com",
	description="Lattice dendrites",
	long_description=long_description,
	long_description_content_type="text/markdown",
	include_package_data=True, # described in MANIFEST.in
	url="https://github.com/urakubo/lattice_dendrites",
	classifiers=[
		"License :: OSI Approved :: MIT License",
		"Programming Language :: Python",
		"Programming Language :: Python :: 3",
		"Programming Language :: Python :: 3.6",
		"Programming Language :: Python :: 3.7",
		"Programming Language :: Python :: 3.8",
		"Programming Language :: Python :: 3.9",
		"Topic :: Scientific/Engineering",
		"Intended Audience :: Science/Research",
		"Operating System :: POSIX",
		"Operating System :: MacOS",
		"Operating System :: Microsoft :: Windows :: Windows 10",
	],
	install_requires=[
		"h5py",
		"numpy>=1.16.1",
		"opencv-python>=4.2.0",
		"pymeshfix>=0.15.0",
		"pyvista>=0.32.1",
		"scikit-image>=0.17.2",
		"trimesh>=3.9.36"
	],
	packages = ['pyLM','pySTDLM','pyLD'],
	package_data={"lm": ['*']},
	python_requires="~=3.6", # >= 3.6 < 4.0
	cmdclass={
	    'develop': PostDevelopCommand,
	    'install': PostInstallCommand,
	},
	entry_points={
		'console_scripts': [
		'check_lm_install_dir = pyLD:check_lm_install_dir',
	],
},

)


