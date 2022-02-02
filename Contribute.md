If you want to contribute to HyperHDG, which is great, we encourage you to re-read the [README](
https://github.com/HyperHDG/HyperHDG/blob/main/README.md) and especially its [Contributions](
https://github.com/HyperHDG/HyperHDG/tree/main#contributions) section. There are several ways
of contributing to HyperHDG:


# Contributing to HyperHDG's code & doxygen: branch `main`

If you want to add contributions to HyperHDG's [main](
https://github.com/HyperHDG/HyperHDG/tree/main) branch, it is important to ensure that its
[push-test](https://github.com/HyperHDG/HyperHDG/blob/main/shell_scripts/push_test.sh) script
succeeds. This is one of two basic requirement for pull-requests to the [main](
https://github.com/HyperHDG/HyperHDG/tree/main) branch to be merged. However, to execute the
script you might need to install `clang-format` (which is used to format all C++ based source code)
and the compilers  specified as `TEST_COMPILER`s in the [push_test.sh](
https://github.com/HyperHDG/HyperHDG/blob/main/shell_scripts/push_test.sh). The other basic
requirement is that the GitHub tests specified in the [workflows](
https://github.com/HyperHDG/HyperHDG/tree/main/.github/workflows) folder succeed. This is
indicated by a green check mark (as opposed to a red cross) on top of the file list of your branch.
A yellow dot indicated that the tests are still running.

With the aforementioned two requirements fulfilled, you can create a pull request which will be
merged to the main branch after an additional revision by one of the [[Authors]] if it can be
automatically merged. You will receive a feedback if any problems with the pull request occurred.

The [doxygen](https://hyperhdg.github.io/auto_pages/doxygen) is automatically generated after
a successful merge to the `main` branch.


# Contributing to HyperHDG's website: branch `gh-pages`

The website of HyperHDG is automatically created from the branch `gh-pages`. Thus, in order to
change the website, please create a pull-request to this branch. Changes will be implemented after
a successful revision by one of the [[Authors]].


# Contributing to HyperHDG's wiki: branch `gh-wiki`

Please do not try to directly change pages within the wiki of HyperHDG. These pages are 
automatically generated from the repository's branch `gh-wiki`. Thus, to change the wiki, please
create a pull-request to this branch. Changes will be implemented after a successful revision by one
of the [[Authors]].


# Other ways of contribution

If you would like to contribute to HyperHDG in a way that is not covered by the aforementioned
points, please contact one of the [[Authors]].
