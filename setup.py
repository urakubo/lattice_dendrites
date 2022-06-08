from setuptools import setup
import os, sys


with open("README.md") as f:
	long_description = f.read()

classifiers = [
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
		"Programming Language :: Python",
		"Topic :: Scientific/Engineering",
		"Intended Audience :: Science/Research"
	]

packages    = ['pyLM','pySTDLM','pyLD','lm']
package_dir = {'pyLM':'pyLM','pySTDLM':'pySTDLM','pyLD':'pyLD'}

major_version = sys.version_info[0] # Major
minor_version = sys.version_info[1] # Minor

print("Python{} on {} is detected.".format(sys.version_info, os.name))
if os.name == 'posix' and major_version == 3 and minor_version == 6:
	classifiers.extend(["Programming Language :: Python :: 3.6",
		"Operating System :: POSIX"])
	package_dir['lm'] = 'lm_py3.6'
	python_requires='3.6'

elif os.name == 'posix' and major_version == 3 and minor_version == 10:
	classifiers.extend(["Programming Language :: Python :: 3.10",
		"Operating System :: POSIX"])
	package_dir['lm'] = 'lm_py3.10'
	python_requires='3.10'

elif major_version == 3:
	classifiers.extend(["Programming Language :: Python :: 3",
		"Operating System :: POSIX",
		"Operating System :: POSIX",
		"Operating System :: Microsoft :: Windows :: Windows 10"])
	package_dir['lm'] = 'lm_null'
	python_requires='>=3.0'

else :
	sys.exit("This version and/or system is not supported.".format(sys.version_info, os.name))


s = setup(
	name="lattice_dendrites",
	version="0.1.0",
	author="Hidetoshi Urakubo",
	author_email="hurakubo@gmail.com",
	description="Lattice dendrites",
	long_description=long_description,
	long_description_content_type="text/markdown",
	include_package_data=True, # described in MANIFEST.in
	url="https://github.com/urakubo/lattice_dendrites",
	classifiers=classifiers,
	install_requires=[
		"h5py",
		"numpy>=1.16.1",
		"opencv-python>=4.2.0",
		"pymeshfix>=0.15.0",
		"pyvista>=0.32.1",
		"scikit-image>=0.17.2",
		"trimesh>=3.9.36"
	],
	packages = packages,
	package_dir = package_dir,
	python_requires=python_requires,
	entry_points={
		'console_scripts': [
		'check_lm_install_dir = pyLD:check_lm_install_dir',
	],
},
)


