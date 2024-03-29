# Welcome to HyperHDG!

[![M2AN Paper](https://img.shields.io/badge/DOI-10.1051%2Fm2an%2F2022011-blue)](
https://doi.org/10.1051/m2an/2022011)

It contains a C++ based library implementing hybrid discontinuous Galerkin methods on extremely
general domains &mdash; that is, all standard (volume) domains, graphs, intersecting surfaces, and
several other types of "domains" that can be interpreted as hypergraphs.

The C++ library can be embedded to C++ programs and/or Python scripts using [Cython](
https://cython.org/). This makes the library a good choice for writing both, high performance codes
and easy to handle teaching scripts/programs.

For more details on HyperHDG, you may visit the [website of HyperHDG](
https://hyperhdg.github.io/HyperHDG), the [Doxygen page](
https://hyperhdg.github.io/auto_pages/doxygen), or its [Wiki](../../wiki). In particular, if you
want to know how to setup HyperHDG, we recommend to visit the [Wiki's setup page](../../wiki/Setup).


## Status of continuous integration

<div align="center">

| Task / Test   | Status                                        | Details                                                   |
|:--------------|:----------------------------------------------|:----------------------------------------------------------|
| Format code   | ![Clang](../../workflows/Clang/badge.svg)     | Check, whether C++ code obeys the formatting rules.       |
| Build library | ![CMake](../../workflows/CMake/badge.svg)     | Build and test library to work for predefined test cases. |
| Make doxygen  | ![Doxygen](../../workflows/Doxygen/badge.svg) | Automatically generate and deploy doxygen to web-page.    |

</div>

The status given in the table refers to the `main` branch. The build test covers building the
library and testing it with several compilers. For more details on the tests, please refer to the
files `.github/workflows`.


# Copyright, License, and Contribution Policy

This directory contains the HyperHDG library.

The HyperHDG library is copyrighted by the HyperHDG authors. This term refers to the people listed
on the [Wiki's page: Authors](../../wiki/Authors).

The HyperHDG library is free software; you can use it, redistribute it, and/or modify it under the 
terms of the <b>GNU Lesser General Public License</b> as published by the Free Software Foundation; 
either <b>version 2.1</b> of the License, or (at your option) any later version. The full text of
the GNU Lesser General Public version 2.1 is quoted in [License.txt](License.txt).


## Contributions

As a contributor to this project, you agree that all of your contributions be governed by the
<b>Developer Certificate of Origin version 1.1</b>. The HyperHDG project does not require copyright
assignments for contributions. This means that the copyright for code contributions in the HyperHDG
project is held by its respective contributors who have each agreed to release their contributed
code under a compatible open source license (LGPL v2.1 for library code). The full text of the 
Developer Certificate of Origin version 1.1 is quoted in [DeveloperCertificateOfOrigin.txt](
DeveloperCertificateOfOrigin.txt).


## Referencing the library

In addition to the terms imposed by the LGPL v2.1 or later, we ask for the following courtesy:

> Every publication presenting numerical results obtained with the help of HyperHDG should state the
> name of the library and cite appropriate HyperHDG references listed at the top of the 
> [Wiki's page: Publications](../../wiki/Publications).

This is the usual, fair way of giving credit to contributors to a scientific result. In addition, it
helps us justify our effort in developing HyperHDG as an academic undertaking.

We keep a list of publications using HyperHDG. Let us know about your publications so that we can 
add them to the aforementioned list. You can do this by emailing the reference information to one of
the principal developers of HyperHDG, cf. [Wiki's authors page](../../wiki/Authors).


## Bundled third party software in the HyperHDG repository

The subdirectory `submodules/` contains third party software. <b>Please note that the software
located there is copyrighted by their respective authors</b> (independent of the HyperHDG authors)
and are covered by different licenses.

The libraries listed above are all open source and as such place few restrictions on their use.
However, if you want to distribute software developed with HyperHDG (in source or binary form) and
you are using the packages above (with source code in `submodules/`), then they may impose
different terms. Please consult the licenses of these packages for more information.

Alternatively, the configuration process of HyperHDG allows you to <i>remove</i> the `submodules/`
directory entirely, or disable any or all of these bundled libraries. In that case, you cannot use
their functionality but you also are not restricted by their license.


## Contact

For further questions regarding licensing and commercial use please contact the principal developers
of HyperHDG, cf. [authors page](../../wiki/Authors), directly.


## Links

- The license can be found in [License.txt](License.txt). It contains the [GNU Lesser General Public
License version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).
- The developer certificate of origin can be found in 
[DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt). It contains the [Developer 
Certificate of Origin version 1.1](https://developercertificate.org/).
- The list of authors and publications can be found on the [authors page](../../wiki/Authors) and
the [publications page](../../wiki/Publications), respectively.
