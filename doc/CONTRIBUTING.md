# CLAIRE: Contributing Guidlines

Go back to [README.md](../README.md).

## Content
* [Overview](#overview)
* [Reporting Bugs and Issues](#bugs)
* [Feature Requests](#features)
* [Contributing to CLAIRE](#contribute)
* [Testing and Benchmarks](#testing)
* [Coding Conventions](#conventions)

## Overview <a name="overview"></a>

We welcome all contributions to CLAIRE. Contributions can be in the form of new code functions, improvements to the documentation, or by pointing out a bug or potential improvement. For general questions or requests, you can contact [Andreas Mang](http://math.uh.edu/~andreas) by email: [andreas [at] math [dot] uh.edu](mailto:andreas@math.uh.edu).

This project and everyone participating in it is governed by the code of conduct found in [doc/CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [Andreas Mang](mailto:andreas@math.uh.edu).



## Reporting Bugs and Issues <a name="bugs"></a>

Bugs and issues can be reported at [https://github.com/andreasmang/claire/issues](https://github.com/andreasmang/claire/issues). Before submitting a ticket, make sure the bug has not already been reported by searching existing issues on GitHub under [issues](https://github.com/andreasmang/claire/issues). We also list some known issues in [doc/README-INSTALL.md](README-INSTALL.md).

If you are reporting a bug or an issue, please include detailed information to help maintainers reproduce the problem.


## Feature Requests <a name="features"></a>

Additional features can be requested at [https://github.com/andreasmang/claire/issues](https://github.com/andreasmang/claire/issues). Before submitting a feature request, make sure the feature has not already been requested by searching existing issues on GitHub under [issues](https://github.com/andreasmang/claire/issues).


## Contributing to CLAIRE <a name="contribute"></a>

To contribute to the software:

1. [Fork](https://docs.github.com/en/free-pro-team@latest/github/getting-started-with-github/fork-a-repo) the repository.
2. Clone the forked repository, add your contributions and push the changes to your fork.
3. Create a [pull request](https://github.com/andreasmang/claire/pulls).


## Testing and Benchmarks <a name="testing"></a>

We have implemented several tests to check the accuracy of our numerical implementation. These are described in more detail in [doc/README-RUNME.md](https://github.com/andreasmang/claire/blob/gpu/doc/README-RUNME.md#testing-and-benchmarks-).


## Coding Conventions <a name="conventions"></a>

Our source code follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html). Please adhere to these coding conventions if you would like to contribute to CLAIRE. Notice that all our routines use the error mechanism implemented in PETSc. We strongly encourage all contributers to include this mechanism to make debugging easier. We have added `tags` to identify the major releases of our software (related mostly to publications for reproducability).
