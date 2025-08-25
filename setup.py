from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        name="effspm._effspm",
        sources=[
            # Python‐binding entrypoint
            "effspm/_effspm.cpp",
            
            "effspm/freq_miner.cpp",
            "effspm/load_inst.cpp",
            "effspm/utility.cpp",
                     # Collector for HTMiner in‐memory patterns
            

            # HTMiner core sources
            
            "effspm/htminer/src/freq_miner.cpp",
            "effspm/htminer/src/load_inst.cpp",
            "effspm/htminer/src/utility.cpp",
            "effspm/htminer/src/build_mdd.cpp",
            
            # BTMiner sources
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
            "effspm",            # for _effspm.cpp
             
            "effspm/htminer/src",       # for HTMiner headers
            "effspm/btminer/src", # for BTMiner headers
            "effspm/largepp/src",
            "effspm/largebm/src",
            "effspm/largehm/src"
        ],
        language="c++",
        extra_compile_args=["-std=c++11", "-O3"],
        
    )
]

setup(
    name="effspm",
    version="0.2.7",
    description="Efficient Sequential Pattern Mining Library",
    author="Yeswanth Vootla",
    packages=["effspm"],
    ext_modules=ext_modules,
    zip_safe=False,
    install_requires=["pybind11"],
)