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
        include_dirs=[pybind11.get_include(), "effspm"],
        language="c++",
        extra_compile_args=["-O3", "-std=c++17"],
    ),
]

setuptools.setup(
    ext_modules=ext_modules,
    zip_safe=False,
)
