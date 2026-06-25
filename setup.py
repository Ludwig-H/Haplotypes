from setuptools import setup, find_packages

setup(
    name="haplo-bench",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pyyaml",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "haplo-bench=haplo_bench.cli:main",
        ],
    },
)
