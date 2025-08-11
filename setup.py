from setuptools import setup, Extension
import platform
import pybind11

is_windows = platform.system() == "Windows"
extra_compile_args = ["/O2", "/std:c++17"] if is_windows else ["-O3", "-std=c++17"]

ext_modules = [
    Extension(
        name="effspm._effspm",
        sources=[
            "effspm/_effspm.cpp",  # binding entrypoint
            "effspm/freq_miner.cpp",
            "effspm/load_inst.cpp",
            "effspm/utility.cpp",

            "effspm/htminer/src/freq_miner.cpp",
            "effspm/htminer/src/load_inst.cpp",
            "effspm/htminer/src/utility.cpp",
            "effspm/htminer/src/build_mdd.cpp",

            "effspm/btminer/src/freq_miner.cpp",
            "effspm/btminer/src/load_inst.cpp",
            "effspm/btminer/src/utility.cpp",
            "effspm/btminer/src/build_mdd.cpp",

            "effspm/largepp/src/freq_miner.cpp",
            "effspm/largepp/src/load_inst.cpp",
            "effspm/largepp/src/utility.cpp",

            "effspm/largebm/src/freq_miner.cpp",
            "effspm/largebm/src/load_inst.cpp",
            "effspm/largebm/src/utility.cpp",
            "effspm/largebm/src/build_mdd.cpp",

            "effspm/largehm/src/freq_miner.cpp",
            "effspm/largehm/src/load_inst.cpp",
            "effspm/largehm/src/utility.cpp",
            "effspm/largehm/src/build_mdd.cpp",
        ],
        include_dirs=[
            pybind11.get_include(),
            "effspm",
            "effspm/htminer/src",
            "effspm/btminer/src",
            "effspm/largepp/src",
            "effspm/largebm/src",
            "effspm/largehm/src",
        ],
        language="c++",
        extra_compile_args=extra_compile_args,
    )
]

setup(
    name="effspm",
    version="0.2.4",
    description="Efficient Sequential Pattern Mining Library",
    author="Yeswanth Vootla",
    packages=["effspm"],
    ext_modules=ext_modules,
    zip_safe=False,
    install_requires=["pybind11>=2.6"],
)
