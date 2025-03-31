# Documentation 🔎
```{rubric} Main Features:
```
- A physics-based model developed in python to estimate rainfall induced by tropical cyclones 🌀.
- Simulate convective tropical cyclone rainfall by correlating the precipitation rate 🌨️ with the total upward wind velocity 🌬️ within the tropical cyclone vortex.
- Integrate with outputs from tropical cyclone downscaling models driven by large-scale environmental conditions from GCMs 🌎.
- Generate detailed spatial-temporal tropical cyclone rainfall patterns that align with hurricane tracks.

This documentation explains how to install and use `pyTCR`.


<!-- Once you've installed, you can use our documentation in two main ways:
1. You jump right in to [notebook examples](examples/examples.md) of how to download TC tracks, simulate TC rainfall, and make visualization.
2. If you prefer to learn about the fundamentals of the library first, you can read about the [mathematical background](math.md). -->

::::{grid} 1 2 2 3
:gutter: 1 1 1 2

:::{grid-item-card} {octicon}`gear;1.5em;sd-mr-1` Installation
:class-body: bg-dark-orange
:link: installation
:link-type: ref

Install *pyTCR* using your favourite package manager.
+++
[Learn more »](./installation.md)
:::

:::{grid-item-card} {octicon}`code-square;1.5em;sd-mr-1` Background
:class-body: bg-dark-green
:link: math
:link-type: ref

Read about the fundamentals of *pyTCR*.
+++
[Learn more »](math.md)
:::

:::{grid-item-card} {octicon}`book;1.5em;sd-mr-1` Tutorials
:class-body: bg-dark-gray
:link: examples/examples
:link-type: doc

Explore jupyter notebooks designed for training *pyTCR*.
+++
[Learn more »](./examples/examples.md)
:::

::::
