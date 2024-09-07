from setuptools import setup, find_packages

setup(
    name='pyTCR',
    version='0.1.0',
    description='Python-based Tropical Cyclone Rainfall',
    author='Phong Le',
    author_email='levuvietphong@gmail.com',
    url='https://github.com/levuvietphong/pyTCR',
    packages=['tcr'],  # Explicitly specify the package to include
    install_requires=[
        'numpy',
        'matplotlib',
        'colorcet',
        'scipy',
        'cartopy',
        'shapely',
        'xarray',
        'netcdf4',
    ],
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)