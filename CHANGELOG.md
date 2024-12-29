#### 1.2 (2024-12-29)

##### Build System / Dependencies

*  update installation method and build packages. - move to pyproject.toml method - build tcr package for import - update notebooks (3c3215c1)

##### Chores

*  restore .gitignore file (ba376ed4)
*  stop tracking .gitignore file (27e2d39d)
*  update .gitignore to exclude chglog template files (c7316996)
*  add changelog capability (36510b13)

##### Documentation Changes

*  update README notebook description (1d39ce9e)
*  update README animation and plots (0643f781)
*  update README and notebook links (02943801)
*  update CHANGELOG for v1.1 (33ee782e)
*  enhance README with example descriptions (dd1f051f)
*  update README file with go-to-top (9c19b4fb)
*  add table of content to README (ab34a92f)
*  update README and fix links to notebooks (a4150544)

##### New Features

*  add functions to load colormaps from cpt files (fd2cc410)
*  switch installation to mamba for speed and efficiency, update README (fbd78441)
*  add color to exdanced plot and rename rm_trks in estimate_radius_wind (ddbdc014)
*  add buffer option to track_landfall (24cec1c4)
*  integrate binder for interactive environment (184f2760)
*  add more examples for demonstrations (34595bac)
*  refine plotting for shapefiles and multiple rain events (b56b0909)
*  enhance plotting with polygon example and exceedance probability (36343196)
*  expand polygon options (22189b97)
*  integrate pyproj and libtiff (cf9668c9)
*  improve data processing and add examples (2121b5e2)
*  refine animations and logo (b7ced2de)
*  include animations and logo (a9039736)
*  enhance GCM I/O functionality (b7ef25d8)
*  add windswathx function and update code format (5e6b6764)

##### Bug Fixes

*  bugs in sfac value in wind.py that increase vertical velocity (a61e1a61)
*  bugs in windprofile to calculate Mm (db691190)
*  a bug in track_landfall that shift longitude 360 deg (76b130ee)
*  handle longitude crossing the prime meridian correctly (9465e8cd)
*  installation bug and move to mamba (f2bc0740)
*  fix bugs in link for notebook example 5 (d4a2a92f)
*  correct model name typo from 2-0 to 1-0 (3572e707)

##### Other Changes

*    - apply autopep8   - eliminate params.py usage   - modify variable suffixes   - include shapefile-derived bbox for extent (dcc79d63)
*    - hurricane tracks from downscaling model   - surface parameter data (bathymetry and drag coeffs)   - load function updates for compatibility (62fd334c)

##### Refactors

*  merge vertical wind functions into one for consistency (230e4ea1)
*  merge vertical wind functions into one for consistency (9c1a9da5)
*  rename functions for enhanced clarity and code comprehension (8ebd8179)
*  optimize codes in terrain_boundary (7b88112f)
*  continue code formatting and update README (31cd1a44)
*  format code with Cursor AI (e5823370)

#### 1.1 (2024-10-08)

##### Build System / Dependencies

*  update installation method and build packages. - move to pyproject.toml method - build tcr package for import - update notebooks (3c3215c1)

##### Chores

*  update .gitignore to exclude chglog template files (c7316996)
*  add changelog capability (36510b13)

##### Documentation Changes

*  enhance README with example descriptions (dd1f051f)
*  update README file with go-to-top (9c19b4fb)
*  add table of content to README (ab34a92f)
*  update README and fix links to notebooks (a4150544)

##### New Features

*  integrate binder for interactive environment (184f2760)
*  add more examples for demonstrations (34595bac)
*  refine plotting for shapefiles and multiple rain events (b56b0909)
*  enhance plotting with polygon example and exceedance probability (36343196)
*  expand polygon options (22189b97)
*  integrate pyproj and libtiff (cf9668c9)
*  improve data processing and add examples (2121b5e2)
*  refine animations and logo (b7ced2de)
*  include animations and logo (a9039736)
*  enhance GCM I/O functionality (b7ef25d8)
*  add windswathx function and update code format (5e6b6764)

##### Bug Fixes

*  fix bugs in link for notebook example 5 (d4a2a92f)
*  correct model name typo from 2-0 to 1-0 (3572e707)

##### Other Changes

*    - apply autopep8   - eliminate params.py usage   - modify variable suffixes   - include shapefile-derived bbox for extent (dcc79d63)
*    - hurricane tracks from downscaling model   - surface parameter data (bathymetry and drag coeffs)   - load function updates for compatibility (62fd334c)

##### Refactors

*  continue code formatting and update README (31cd1a44)
*  format code with Cursor AI (e5823370)

