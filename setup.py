import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lattice_dendrites",
    version="0.1.0",
    author="hidetoshi-urakubo",
    author_email="hurakubo@gmail.com",
    description="Lattice dendrite",
    long_description=read('README.md'),
    long_description_content_type="text/markdown",
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
		"json",
		"numpy>=1.16.1",
		"opencv-python>=4.5.4.60",
		"pymeshfix>=0.15.0",
		"pyvista>=0.32.1",
		"scikit-image>=0.17.2",
		"trimesh>=3.9.36"
	],
    packages = ['pyLM','pySTDLM','pyLD'],
    entry_points = {
        'console_scripts': ['sample_command = sample_command.sample_command:main']
    },
    python_requires="~=3.6", # >= 3.6 < 4.0
)
