# Change Log

## [pyTCR-1.2.2](https://github.com/levuvietphong/pyTCR/compare/pyTCR-1.2.1..pyTCR-1.2.2) - 2025-04-09

### 🐛 Bug Fixes

- Resolve a bug in iodata to get bounding box - ([765cb5d](https://github.com/levuvietphong/pyTCR/commit/765cb5dca6a044ae4cafdac24bc990431a6edcad))

### 📚 Documentation

- Update README file and change docs/paper to avoid confusion for paper compilation - ([1ee90b1](https://github.com/levuvietphong/pyTCR/commit/1ee90b134e7a5d8c59264b706c18bd3c654dd261))
- Improve API docstring for I/O and plotting functions - ([c7d9430](https://github.com/levuvietphong/pyTCR/commit/c7d9430bff8b19fe3e99aed3038ae81698ffc09e))
- Add conda installation to README - ([c9b15a4](https://github.com/levuvietphong/pyTCR/commit/c9b15a4415b7e01f2fa78a41dd010381dd9422f5))
- Add conda installation to README - ([13b7ac8](https://github.com/levuvietphong/pyTCR/commit/13b7ac82e63c7bcbb1620e02c8433d9cb03e604b))
- Update version 1.2.2 development and CI - ([016a043](https://github.com/levuvietphong/pyTCR/commit/016a04367dddcb0d19b981661011e5bad55a8849))

### ⚙️ Miscellaneous Tasks

- Add workflow to check the paper in joss - ([ff98a19](https://github.com/levuvietphong/pyTCR/commit/ff98a195d4aaead447f9f1169e829847819f139f))

## [pyTCR-1.2.1](https://github.com/levuvietphong/pyTCR/compare/pyTCR-1.2..pyTCR-1.2.1) - 2025-04-03

### ⛰️  Features

- Migrate installation to PyPI and add contribution guides - ([1742f22](https://github.com/levuvietphong/pyTCR/commit/1742f222dbd019ec725187e641b0ad73ceafce3b))
- Add sphinx read-the-docs for documentation - ([96813b0](https://github.com/levuvietphong/pyTCR/commit/96813b0c1c9c9dd9b66516c80061150fbed978c3))
- Add config and requirements for read-the-docs - ([bf4fe90](https://github.com/levuvietphong/pyTCR/commit/bf4fe909fd07bdb1e3d3223b33220c12a0db540d))
- Add documentation using jupyter book - ([3e8e568](https://github.com/levuvietphong/pyTCR/commit/3e8e568db1c3696cc5495ffa1ab18ceedfe6e428))
- Add pytest for various tcr functions - ([3d47e5c](https://github.com/levuvietphong/pyTCR/commit/3d47e5c1037832071429318fa424be503ce14b67))

### 🐛 Bug Fixes

- Update units in all functions and docs - ([a0b1155](https://github.com/levuvietphong/pyTCR/commit/a0b11558d0bd77387a705270b50496919c3e17ac))
- Update version in pyproject - ([d120d77](https://github.com/levuvietphong/pyTCR/commit/d120d772ec0389c2ef75453ee0aff7fa8cb8793a))
- Resolve indexing issue with short velocity time series - ([5f3b51a](https://github.com/levuvietphong/pyTCR/commit/5f3b51a965661994732c74cb2d4c2ea5ca4dd035))
- Modify CI for PyPI migration - ([d7b3c6c](https://github.com/levuvietphong/pyTCR/commit/d7b3c6c8f261f723bc21b3185e59cb37a0d1b217))

### 📚 Documentation

- Update documentation website and notebooks - ([69617c1](https://github.com/levuvietphong/pyTCR/commit/69617c1a00cad0493b93c62308dcb53999e6530a))
- Update readthedocs yaml file for sphinx - ([4d4ec33](https://github.com/levuvietphong/pyTCR/commit/4d4ec33c845bcafa5bc4b8eac1af2fd22b049009))
- Add README for PyPI description - ([fc588ba](https://github.com/levuvietphong/pyTCR/commit/fc588bab2ba652ac3b9c4dc0368dcca9e8303190))
- Fix bibliography and disable external execution - ([645bf93](https://github.com/levuvietphong/pyTCR/commit/645bf93d9218a3e65a23148a24a0bb5de206145f))
- Fix bibliography and disable external execution - ([5cf7bb2](https://github.com/levuvietphong/pyTCR/commit/5cf7bb26fe22b55e52a220b2e058a07abdfb9fc3))
- Update mathematical background and references - ([f0e0474](https://github.com/levuvietphong/pyTCR/commit/f0e047423fe46d8abec9e8929ee7cc709bae2099))
- Revise documentation website - ([f2ac246](https://github.com/levuvietphong/pyTCR/commit/f2ac24644a7ab4e19db4518adab0d1c6e2665376))
- Update docstring and api function of dirdata - ([e083d1a](https://github.com/levuvietphong/pyTCR/commit/e083d1a25a80ce7ea86c42bf43d5551a4b5f64fc))
- Minor update on the introduction and version - ([62b2d3e](https://github.com/levuvietphong/pyTCR/commit/62b2d3e9e2a863e793d98e86908ce7b33d716456))
- Move changelog to git-cliff - ([5c10c00](https://github.com/levuvietphong/pyTCR/commit/5c10c00a90cd9b81013272e06ed86f1f4106b37c))

### ⚙️ Miscellaneous Tasks

- Remove matlab test files - ([ea279c1](https://github.com/levuvietphong/pyTCR/commit/ea279c14d229200af1da888f8660da6bccef01a0))

## New Contributors ❤️

* @levuvietphong made their first contribution
* @ecoon made their first contribution in [#6](https://github.com/levuvietphong/pyTCR/pull/6)
## [pyTCR-1.2](https://github.com/levuvietphong/pyTCR/compare/pyTCR-1.1..pyTCR-1.2) - 2024-12-29

### ⛰️  Features

- Add buffer option to track_landfall - ([24cec1c](https://github.com/levuvietphong/pyTCR/commit/24cec1c4073366dcf96d074dfff96e2c4341ba21))
- Add color to exdanced plot and rename rm_trks in estimate_radius_wind - ([ddbdc01](https://github.com/levuvietphong/pyTCR/commit/ddbdc01452f876762229174d450989921e0f85b8))
- Switch installation to mamba for speed and efficiency, update README - ([fbd7844](https://github.com/levuvietphong/pyTCR/commit/fbd784416ba6f06df3283431cda6ee6de43678e8))
- Add functions to load colormaps from cpt files - ([fd2cc41](https://github.com/levuvietphong/pyTCR/commit/fd2cc410c2fa6cf875aa4937f0494bd0dadf8f78))

### 🐛 Bug Fixes

- Handle longitude crossing the prime meridian correctly - ([9465e8c](https://github.com/levuvietphong/pyTCR/commit/9465e8cd31f111a1a943db739a64e7637f2e0a57))
- A bug in track_landfall that shift longitude 360 deg - ([76b130e](https://github.com/levuvietphong/pyTCR/commit/76b130ee485b4bc9648f9159ea57605c7aab43d1))
- Installation bug and move to mamba - ([f2bc074](https://github.com/levuvietphong/pyTCR/commit/f2bc0740359b73988ce905c7383b681c7e3e5469))
- Bugs in windprofile to calculate Mm - ([db69119](https://github.com/levuvietphong/pyTCR/commit/db691190b8c7d881b3695d24b5007feaacbf6191))
- Bugs in sfac value in wind.py that increase vertical velocity - ([a61e1a6](https://github.com/levuvietphong/pyTCR/commit/a61e1a615a053d5de7ef06adcbf9ae4aeb30b339))

### 🚜 Refactor

- Optimize codes in terrain_boundary - ([7b88112](https://github.com/levuvietphong/pyTCR/commit/7b88112fa0d6035e6b6e5a88b6845e6b3d2f2c14))
- Rename functions for enhanced clarity and code comprehension - ([8ebd817](https://github.com/levuvietphong/pyTCR/commit/8ebd8179fe7ec95bcf66a24888e5de90cfa162a7))
- Merge vertical wind functions into one for consistency - ([9c1a9da](https://github.com/levuvietphong/pyTCR/commit/9c1a9da59cdbfcb2efaa8fd13a6327150d1bcecf))
- Merge vertical wind functions into one for consistency - ([230e4ea](https://github.com/levuvietphong/pyTCR/commit/230e4ea18a47d605655d439be7c568a6114055ce))

### 📚 Documentation

- Update README and notebook links - ([0294380](https://github.com/levuvietphong/pyTCR/commit/02943801ee34d6113bfd612706e7285b7ffdea6a))
- Update README animation and plots - ([0643f78](https://github.com/levuvietphong/pyTCR/commit/0643f7819280b60ad8dc63a0990d3e986fac2812))
- Update README notebook description - ([1d39ce9](https://github.com/levuvietphong/pyTCR/commit/1d39ce9ed39b6f38dba3233e8b85b7c06fd5b9e1))
- Update CHANGELOG.md for ver1.2 - ([a31d687](https://github.com/levuvietphong/pyTCR/commit/a31d6876ee049e03268120bdf4adc2f9506db86c))

### ⚙️ Miscellaneous Tasks

- Stop tracking .gitignore file - ([27e2d39](https://github.com/levuvietphong/pyTCR/commit/27e2d39d73366d255e0cef17cc6f7bb895c2ee6b))
- Restore .gitignore file - ([ba376ed](https://github.com/levuvietphong/pyTCR/commit/ba376ed4bb4914c6c3aa971fafe67491a009f0ea))

## [pyTCR-1.1](https://github.com/levuvietphong/pyTCR/compare/pyTCR-1.0..pyTCR-1.1) - 2024-10-08

### ⛰️  Features

- Add windswathx function and update code format - ([5e6b676](https://github.com/levuvietphong/pyTCR/commit/5e6b67640dd80347bf44c92e771f156e755cad3e))
- Enhance GCM I/O functionality - ([b7ef25d](https://github.com/levuvietphong/pyTCR/commit/b7ef25d83eba81287cefe22213d90a1e6e66f96e))
- Include animations and logo - ([a903973](https://github.com/levuvietphong/pyTCR/commit/a9039736105ea5f75a685edffd0d92f8a8a21738))
- Refine animations and logo - ([b7ced2d](https://github.com/levuvietphong/pyTCR/commit/b7ced2deb1e11e4bf74bb94a26b2f1a79a4cf2cd))
- Improve data processing and add examples - ([2121b5e](https://github.com/levuvietphong/pyTCR/commit/2121b5e2f31dc8e09d258310e2121eb1033e84ff))
- Integrate pyproj and libtiff - ([cf9668c](https://github.com/levuvietphong/pyTCR/commit/cf9668c99aa5e8aef8bbbb4d7292a123437747f5))
- Expand polygon options - ([22189b9](https://github.com/levuvietphong/pyTCR/commit/22189b97a07f8bc51c792d7cdb825b955d598aeb))
- Enhance plotting with polygon example and exceedance probability - ([3634319](https://github.com/levuvietphong/pyTCR/commit/363431960093d8d2c82c10864eb24fe2ea6a96e6))
- Refine plotting for shapefiles and multiple rain events - ([b56b090](https://github.com/levuvietphong/pyTCR/commit/b56b0909f47c717ece2e2a753a3fe859ca3d59fc))
- Add more examples for demonstrations - ([34595ba](https://github.com/levuvietphong/pyTCR/commit/34595bac6dc97b8edad4f79b25091b7abbc58b56))
- Integrate binder for interactive environment - ([184f276](https://github.com/levuvietphong/pyTCR/commit/184f2760905b1adea72a5452d04cd015acb9c868))

### 🐛 Bug Fixes

- Correct model name typo from 2-0 to 1-0 - ([3572e70](https://github.com/levuvietphong/pyTCR/commit/3572e7077604961bb7de7f96a31016d631249e6b))
- Fix bugs in link for notebook example 5 - ([d4a2a92](https://github.com/levuvietphong/pyTCR/commit/d4a2a92fe9892f734e208b43a20bff9f72864649))

### 🚜 Refactor

- Format code with Cursor AI - ([e582337](https://github.com/levuvietphong/pyTCR/commit/e5823370e022f2656178f91e7a221220590a9cf5))
- Continue code formatting and update README - ([31cd1a4](https://github.com/levuvietphong/pyTCR/commit/31cd1a4456810596231d72c630ba8c9e13969793))

### 📚 Documentation

- Update README and fix links to notebooks - ([a415054](https://github.com/levuvietphong/pyTCR/commit/a41505448001943ada8bf788e36a2214547fbb3f))
- Add table of content to README - ([ab34a92](https://github.com/levuvietphong/pyTCR/commit/ab34a92f1e0b25ea656d3cfff5e011290faf7413))
- Update README file with go-to-top - ([9c19b4f](https://github.com/levuvietphong/pyTCR/commit/9c19b4fb525f036ecf94ddc88c05eb47c003eded))
- Enhance README with example descriptions - ([dd1f051](https://github.com/levuvietphong/pyTCR/commit/dd1f051fcffaba7ccc8d579e556cd7a6202f9185))
- Update CHANGELOG for v1.1 - ([33ee782](https://github.com/levuvietphong/pyTCR/commit/33ee782e54062929c585895fc5cffb5e090e62b9))

### ⚙️ Miscellaneous Tasks

- Add changelog capability - ([36510b1](https://github.com/levuvietphong/pyTCR/commit/36510b1395745294737aa87b841d3ae33a04c880))
- Update .gitignore to exclude chglog template files - ([c731699](https://github.com/levuvietphong/pyTCR/commit/c73169967238c1a4b4985fbf61d250cc1584a522))

<!-- generated by git-cliff -->
