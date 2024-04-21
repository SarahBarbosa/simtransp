from setuptools import setup

setup(
    name="simtransp", 
    version = "1.0.0",
    description = "Simulador de trânsito planetário simples",
    long_description = "README.md",
    author = "Sarah G. A. Barbosa",
    author_email = "sarahg.aroucha@gmail.com",
    packages = ["simtransp"],
    install_requires = ["numpy", "PyAstronomy", "astropy", "matplotlib", "numpy"], 
    zip_safe = False,
    )