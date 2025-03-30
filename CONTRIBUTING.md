<!-- omit in toc -->
# pyTCR Contribution Guide

Thank you for your interest in contributing to pyTCR. Contributions are greatly appreciated!

### Getting Started
Ready to contribute? Here's how to set up *pyTCR* for local development.

1. Fork the [pyTCR](https://github.com/levuvietphong/pyTCR) repo on GitHub.

2. Clone your fork locally:
    ```sh
    git clone https://github.com/your-username/pyTCR.git
    ```
3. We recommend using `mamba` to manage dependencies. If you already have `mamba` installed, skip to the next step. If not, you can install it by running:
    ```sh    
    conda install -n base -c conda-forge mamba
    ```
4. Install your local copy into a conda environment:
    ```sh
    mamba env create -f environment.yml
    conda activate pyTCR
    pip install -e .
    ```

4. Create a branch for local development:
    ```sh
    git checkout -b feature/your-feature-name
    ```
    Now you can make your changes locally.


5. Commit your changes and push your branch to GitHub:
    ```sh
    git add .
    git commit -m "Add feature: Your detailed feature description"
    git push origin feature/your-feature-name
    ```
6. Submit a pull request through the GitHub website.

### Report Bugs
Report bugs at https://github.com/levuvietphong/pyTCR/issues.

If you are reporting a bug, please include:

- Your operating system name and version.
- Any details about your local setup that might be helpful in troubleshooting.
- Detailed steps to reproduce the bug.

### Pull Request Guidelines
Before you submit a pull request, check that it meets these guidelines:

- The pull request should include tests.
- If the pull request adds functionality, the docs should be updated. Put your new functionality into a function with a docstring, and add the feature to the list in README.

### Code Style
This project follows PEP 8 coding conventions.

### License
By contributing, you agree that your contributions will be licensed under the same license as this project. Thank you for your contributions! 🚀
