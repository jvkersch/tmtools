import glob
from setuptools import setup, find_packages

from pybind11.setup_helpers import Pybind11Extension


ext_modules = [
    Pybind11Extension(
        "tmtools._bindings",
        sorted(glob.glob("src/*.cpp")) + ["src/extern/TMalign-modified.cpp"],
        cxx_std=14,
    ),
]

setup(
    name="tmtools",
    version="0.0.1",
    author="Joris Vankerschaver",
    author_email="joris.vankerschaver@gmail.com",
    url="https://github.com/jvkersch/tmtools",
    description="Python bindings around the TM-align code",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
    ],
    license="GPLv3",
    platforms=["Linux", "Mac OS-X", "Unix"],
    ext_modules=ext_modules,
    packages=find_packages(),
    package_data={
        "tmtools": ["data/*"],
    },
    include_package_data=True,
    install_requires=[
        "biopython",
        "numpy",
    ],
)
