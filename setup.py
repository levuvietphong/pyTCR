from setuptools import setup, find_packages

setup(
    name='pyTCR',  # Replace with your package name
    version='0.1.0',  # Initial version
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'colorcet',
        'scipy',
        'cartopy',
        'shapely',
        'xarray',
        'netcdf4',
        'xarray',
        'shapely',
    ],
    include_package_data=True,
    description='Python-based Tropical Cyclone Rainfall',
    author='Phong Le',
    author_email='levuvietphong@gmail.com',
    url='https://github.com/levuvietphong/pyTCR',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
