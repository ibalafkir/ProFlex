from setuptools import setup, find_packages

setup(
    name="proflex",
    author="Ismael Balafkir, Anna M. Diaz-Rovira",
    author_email="ibalafkir@gmail.com, annadiarov@gmail.com",
    description=("Explore ProFlex protocol for rigid-body protein-protein "
                 "docking, interface diffusion, and Monte Carlo rotations/"
                 "translations using PyDock, RFDiffusion, and PELE."),
    # license="MIT",
    keywords="protein-protein interaction, protein-protein docking, "
             "interface diffusion, Monte Carlo, PyDock, RFDiffusion, PELE",
    url="https://github.com/ibalafkir/ProFlex",
    packages=find_packages(exclude=('tests', 'docs')),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)