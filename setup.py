import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lattice_dendrites",
    version="0.1.0",
    author="hidetoshi-urakubo",
    author_email="hurakubo@gmail.com",
    description="You can receive the message 'Hello!!!'",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/urakubo/lattice_dendrites",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages = ['pyLM','pySTDLM','pyLD'],
    entry_points = {
        'console_scripts': ['sample_command = sample_command.sample_command:main']
    },
    python_requires='>=3.6',
)
