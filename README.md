# What is HyperHDG?

This repository contains the code to write a hybrid discontinuous Galerkin based solver for PDEs
defined on hypergraphs.


## How to start?

Before you start using HyperHDG, you need to install some packages. Having Ubuntu 20.04 LTS as
operating system this can be done using the commands:

```
$ sudo apt-get update
$ sudo apt-get dist-upgrade
$ sudo apt-get autoremove
$ sudo apt-get install git python3-dev python3-numpy python3-scipy cython libblas-dev liblapack-dev
```

The first three lines bring your operating system into a clear state. The last line installs the
packages that are necessary to fully use HyperHDG. Afterwards, you can start obtaining HyperHDG.

To do so, enter the directory you want to clone HyperHDG into and do the following steps:

1. Clone the repository and give it the name "your_name":

   1. To clone this repository with https use
      
          $ git clone https://github.com/AndreasRupp/HyperHDG.git your_name

   2. To clone this repository with ssh use
      
          $ git git@github.com:AndreasRupp/HyperHDG.git your_name

2. Enter directory using

       $ cd your_name

3. Execute the script `setup.sh` to install HyperHDG by

       $ chmod 700 setup.sh
       $ ./setup.sh



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
Developer Certificate of Origin version 1.1 is quoted in [DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt).


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