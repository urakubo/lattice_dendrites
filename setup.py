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
		print(self.install_lib)

		filename    = os.path.join(self.install_lib, 'lm', "lm_py3.6_cuda11.0.tar.gz")
		extract_dir = os.path.join(self.install_lib, 'lm')
		with tarfile.open(filename, "r:*") as f:
			f.extractall(path=extract_dir, numeric_owner=True)
		major_version = sys.version_info[0] # Major
		minor_version = sys.version_info[1] # Minor

		print('Processed here.')
		if os.name== 'posix' and major_version == 3 and minor_version == 6:
			cuda_version = locate_cuda()
			if cuda_version == '11.0':
				print('Postprocess: posix (linux, mac), python3.6, cuda11.0.')
				print('Postprocess: lm_py3.6_cuda11.0 is extracted.')
				target_file = os.path.join(self.install_lib, 'lm', "lm_py3.6_cuda11.0.tar.gz")
				extract_dir = os.path.join(self.install_lib, 'lm')
				with tarfile.open(target_file, "r:*") as f:
					f.extractall(path=extract_dir, numeric_owner=True)
				print('Please set the following two paths.')
				print('Please export PATH={}:$PATH',  os.path.join(extract_dir,'bin') )
				print('Please export LD_LIBRARY_PATH={}:$LD_LIBRARY_PATH',  os.path.join(extract_dir,'lib') )


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
	description="Lattice dendrite",
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
	packages = ['pyLM','pySTDLM','pyLD','lm'],
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


