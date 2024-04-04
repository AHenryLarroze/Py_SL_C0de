import setuptools

VERSION = "1.0.0"  # PEP-440

NAME = "slcode"

INSTALL_REQUIRES = [
    "joblib==1.2.0",
    "matplotlib==3.7.0",
    "numpy==1.24.2",
    "plotly==5.13.0",
    "pyshtools==4.10.1",
    "scipy==1.10.0",
    "sphinx_rtd_theme==1.2.0",
    "tomli==2.0.1",
    "stripy==2.2",
    "netcdf4==1.6.2",
    "jupyter==1.0.0",
]


setuptools.setup(
    name=NAME,
    version=VERSION,
    description="SL_C0de is a python library based on the matlab code SL_code. This module is used to calculate GIA over the world but alse the sediment isostasy. It include function to calculate radial symetric auto gravitational earth visco-elastic response. This module also include the possibility to resolve the Sea level equation based on the publication of Mitrovica and Austermann ... ajouter des refs. The aim of this library is to be easy to include new sources of load and link to other python library.",
    url="https://github.com/AHenryLarroze/Py_SL_C0de",
    project_urls={
        "Homepage": "https://py-sl-c0de.readthedocs.io/en/latest/",
        "Source Code": "https://github.com/AHenryLarroze/Py_SL_C0de",
    },
    author="Adrien Henry-Larroze",
    author_email="adrien_henry@orange.fr",
    license="GPL-3.0",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    # Snowpark requires Python 3.8
    python_requires=">=3.7",
    # Requirements
    install_requires=INSTALL_REQUIRES,
    packages=["slcode"],
    long_description=open("README.rst").read(),
    long_description_content_type="text/markdown",
)