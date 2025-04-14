# setup.py
import setuptools
from setuptools import Extension
import pybind11

ext_modules = [
    Extension(
        "effspm._core",
        sources=[
            "effspm/_core.cpp",
            "effspm/load_inst.cpp",
            "effspm/freq_miner.cpp",
            "effspm/utility.cpp",
        ],
        include_dirs=[
            pybind11.get_include(),
            "effspm",
        ],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17"],
    ),
]

setuptools.setup(
    name="effspm",
    version="0.1.0",  # keep in sync with pyproject.toml
    author="yeshu999",
    author_email="vootlayeswanth20@gmail.com",
    description="Prefix‑Projection sequential pattern mining",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yeshu999/effspm",
    packages=setuptools.find_packages(where="."),
    ext_modules=ext_modules,
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=["pybind11>=2.6"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
