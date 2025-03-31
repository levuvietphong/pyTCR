# PyTCR: A tropical cyclone rainfall model for python

![](https://img.shields.io/github/license/levuvietphong/pyTCR)
![](https://img.shields.io/github/issues/levuvietphong/pyTCR)
![](https://img.shields.io/github/forks/levuvietphong/pyTCR)
![](https://img.shields.io/github/last-commit/levuvietphong/pyTCR)
![](https://img.shields.io/github/v/release/levuvietphong/pyTCR)
![](https://img.shields.io/github/actions/workflow/status/levuvietphong/pyTCR/CI.yml)

pyTCR is a physics-based model developed in python to estimate rainfall induced by tropical cyclones (TCs). It is largely based on the TCR model described by [Zhu *et al.*, 2013](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013GL058284) and [Lu *et al.*, 2018](https://journals.ametsoc.org/view/journals/atsc/75/7/jas-d-17-0264.1.xml). PyTCR simulates convective TC rainfall by correlating the precipitation rate with the total upward velocity within the TC vortex. It integrates seamlessly with outputs from [a tropical cyclone downscaling model](https://github.com/linjonathan/tropical_cyclone_risk) (see [Lin *et al.,* 2023](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023MS003686)) to generate detailed spatial-temporal rainfall patterns that align with hurricane tracks.

## Installation

`pyTCR` can be installed from `PyPI` using:

```sh
pip install pyTCR
```

## Sources
The source code is available on [GitHub](https://github.com/levuvietphong/pyTCR) under the [MIT Licence](https://github.com/levuvietphong/pyTCR/blob/main/LICENSE).
