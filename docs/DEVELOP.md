# Development guide

## Releasing a new version

1. Integrate all required changes on the `master` branch.
2. Increment the version in [setup.py](../setup.py) on the `master` branch.
3. Ensure that the tests pass: [![.github/workflows/build-and-test.yaml](https://github.com/quadprog/quadprog/actions/workflows/build-and-test.yaml/badge.svg?branch=master)](https://github.com/quadprog/quadprog/actions/workflows/build-and-test.yaml).
4. Download the built artifact from the workflow. The artifact is named `sdist` and the download will be a zipfile containing a .tar.gz.
5. Upload the .tar.gz to PyPI.
6. Publish a new Release in the GitHub repo, creating a new tag, and including the PyPI release in the description.