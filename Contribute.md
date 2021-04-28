If you want to contribute to HyperHDG, which is great, we encourage you to re-read the [README](
https://github.com/AndreasRupp/HyperHDG/blob/master/README) and especially its [Contributions](
https://github.com/AndreasRupp/HyperHDG/tree/master#contributions) section. There are several ways
of contributing to HyperHDG:


# Contributing to HyperHDG's code & doxygen: branch `master`

If you want to add contributions to HyperHDG's [master](
https://github.com/AndreasRupp/HyperHDG/tree/master) branch, it is important to ensure that its
[push-test](https://github.com/AndreasRupp/HyperHDG/blob/master/shell_scripts/push_test.sh) script
succeeds. This is one of two basic requirement for pull-requests to the [master](
https://github.com/AndreasRupp/HyperHDG/tree/master) branch to be merged. However, to execute the
script you might need to install `clang-format` (which is used to format all C++ based source code)
and the compilers  specified as `TEST_COMPILER`s in the [Makefile](
https://github.com/AndreasRupp/HyperHDG/blob/master/Makefile). The other basic requirement is that
the GitHub tests specified in the [workflows](
https://github.com/AndreasRupp/HyperHDG/tree/master/.github/workflows) folder succeed. This is
indicated by a green check mark (as opposed to a red cross) on top of the file list of your branch.
A yellow dot indicated that the tests are still running.

With the aforementioned two requirements fulfilled, you can create a pull request which will be
merged to the master branch after an additional revision by one of the [Authors](Authors) if it can
be automatically merged. You will receive a feedback if any problems with the pull request occurred.

The [doxygen](https://hyperhdg.github.io/auto_pages/doxygen) is automatically generated after
a successful merge to the `master` branch.


# Contributing to HyperHDG's website: branch `gh-pages`

The website of HyperHDG is automatically created from the branch `gh-pages`. Thus, in order to
change the website, please create a pull-request to this branch. Changes will be implemented after
a successful revision by one of the [Authors](Authors).


# Contributing to HyperHDG's wiki: branch `gh-wiki`

Please do not try to directly change pages within the wiki of HyperHDG. These pages are 
automatically generated from the repository's branch `gh-wiki`. Thus, to change the wiki, please
create a pull-request to this branch. Changes will be implemented after a successful revision by one
of the [Authors](Authors).


# Other ways of contribution

If you would like to contribute to HyperHDG in a way that is not covered by the aforementioned
points, please contact one of the [Authors](Authors).
