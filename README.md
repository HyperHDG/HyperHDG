# What is HyperHDG?

This repository contains the code to write a hybrid discontinuous Galerkin based solver for PDEs
defined on hypergraphs. 


# Welcome to HyperHDG!

It contains a C++ based library implementing hybrid discontinuous Galerkin methods on extremely
general domains &mdsh; that is, all standard (volume) domains, graphs, self-intersecting surfaces,
and several other types of "domains" that can be interpreted as hypergraphs.

The C++ library can be embedded to C++ programs and/or Python scripts. This makes the library a good
choice for writing both, high performance codes and easy to handle teaching scripts / programs.

For more details on HyperHDG, you may visit

- the [website of HyperHDG](https://andreasrupp.github.io/HyperHDG),
- the [Doxygen page](https://andreasrupp.github.io/HyperHDG_pages/doxygen),
- or its [Wiki](https://github.com/AndreasRupp/HyperHDG/wiki) containing guides on how to [setup](
https://github.com/AndreasRupp/HyperHDG/wiki/Setup) and [use](
https://github.com/AndreasRupp/HyperHDG/wiki/Usage) the library.


# Copyright, License, and Contribution Policy

This directory contains the HyperHDG library.

The HyperHDG library is copyrighted by the HyperHDG authors. This term refers to the people listed
in the file [Authors.txt](Authors.txt).

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
> name of the library and cite appropriate HyperHDG references listed at the top of the file
> [Publications.txt](Publications.txt).

This is the usual, fair way of giving credit to contributors to a scientific result. In addition, it
helps us justify our effort in developing HyperHDG as an academic undertaking.

We keep a list of publications using HyperHDG. Let us know about your publications so that we can 
add them to the aforementioned list. You can do this by emailing the reference information to one of
the principal developers of HyperHDG, cf. [Authors.txt](Authors.txt).


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
of HyperHDG, cf. [Authors.txt](Authors.txt), directly.


## Links

- The license can be found in [License.txt](License.txt). It contains the [GNU Lesser General Public
License version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).
- The developer certificate of origin can be found in 
[DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt). It contains the [Developer 
Certificate of Origin version 1.1](https://developercertificate.org/).
- The list of authors and publications can be found in [Authors.txt](Authors.txt) and 
[Publications.txt](Publications.txt), respectively.
