(installation)=

# Installation

To install `pyTCR`, choose your preferred installation method:
````{tab} PyPI
To install using `pip`, run:
```bash
$ pip install pyTCR
```
````

````{tab} Conda Forge
To install using `conda`, run:
```bash
$ conda install -c conda-forge pyTCR
```
````

````{tab} Source code

Ensure that you have [`mamba`](https://github.com/mamba-org/mamba) installed.
1. Clone the repository:
```bash
$ git clone https://github.com/levuvietphong/pyTCR.git
$ mamba env create -f environment.yml
$ conda activate pyTCR
```

2. Build from source:
```sh
$ cd pyTCR
$ pip install .
```
````

```{note}
Python 3.11 or higher is required to install `pyTCR`. To verify that the installation was successful, check the `pyTCR` version with:

    $ python -c "import tcr; print(tcr.__version__)"
    
If you have any difficulties during installation, please [open an issue here](https://github.com/levuvietphong/pyTCR/issues). <br>
If you want to contribute to `pyTCR`, see the [contributing guide](https://github.com/levuvietphong/pyTCR/blob/main/CONTRIBUTING.md).

```
