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


def _read_long_description():
    with open("README.md", encoding="utf-8") as fp:
        return fp.read()


setup(
    name="tmtools",
    version="0.2.0",
    author="Joris Vankerschaver",
    author_email="joris.vankerschaver@gmail.com",
    url="https://github.com/jvkersch/tmtools",
    long_description=_read_long_description(),
    long_description_content_type='text/markdown',
    description="Python bindings around the TM-align code for structural alignment of proteins",  # noqa
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering",
    ],
    license="GPLv3",
    platforms=["Linux", "Mac OS-X", "Unix", "Windows"],
    ext_modules=ext_modules,
    packages=find_packages(),
    package_data={
        "tmtools": ["data/*"],
    },
    include_package_data=True,
    install_requires=[
        "numpy",
    ],
)
