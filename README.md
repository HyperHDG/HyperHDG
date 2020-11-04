# What is HyperHDG?

This repository contains the code to write a hybrid discontinuous Galerkin based solver for PDEs
defined on hypergraphs.


## How to start?

To start with this repository enter the directory you want to clone this repository into and do the
following steps:

1. Clone the repository and give it the name "your_name":

   1. To clone this repository with https use:
      
          $ git clone https://github.com/AndreasRupp/HyperHDG.git your_name

   2. To clone this repository with ssh use:
      
          $ git git@github.com:AndreasRupp/HyperHDG.git your_name

2. Enter directory using:

       $ cd your_name

3. Initialize submodules using:

       $ git submodule update --init --recursive



# HyperHDG Copyright, License, and Developer Certificate of Origin

This directory contains the hyperHDG library.

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
Developer Certificate of Origin version 1.1 is quoted in [DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt).


## Referencing the library

In addition to the terms imposed by the LGPL v2.1 or later, we ask for the following courtesy:

> Every publication presenting numerical results obtained with the help of hyperHDG should state the
> name of the library and cite appropriate hyperHDG references listed at the top of the file
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
However, if you want to distribute software developed with deal.II (in source or binary form) and
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
- The developer certificate of origign can be found in 
[DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt). It contains the [Developer 
Certificate of Origin version 1.1](https://developercertificate.org/).
- The list of authors can be found in [Authors.txt](Authors.txt).