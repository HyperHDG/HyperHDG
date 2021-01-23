PROJECT     	= HyperHDG
.PHONY:       	build clean distclean clean_build clean_domains clean_doxygen clean_output \
                clean_pycache doxygen format submodules test_all_compilers test_compiler

# List of compilers that HyperHDG is tested to run with.
TEST_COMPILER = clang++-10 clang++-11 g++-10


####################################################################################################
# The default compiler of GNUmake is CXX=g++. This can be changed by
# an environment variable of the same name or by the command line
# argument CXX=compiler when running make
#
# The respective compiler must support C++20. If $(CXX) is listed in
# $(TEST_COMPILER), the code is pretty safe to run, since it is tested against the usage of all
# compilers defined as TEST_COMPILER.
####################################################################################################

####################################################################################################
# The default target generates a build subdirectory and then configures the library there, builds
# it and runs the tests
####################################################################################################

build:
	mkdir -p build
	cd build; cmake -DCMAKE_CXX_COMPILER=$(CXX) ..
	cd build; make
	cd build; make test


## Clean up the automatically generated files of the library.
clean:
	$(MAKE) clean_build
	$(MAKE) clean_domains
	$(MAKE) clean_doxygen
	$(MAKE) clean_pycache


## Clean up the automatically generated files of the library and output files.
distclean:
	$(MAKE) clean
	$(MAKE) clean_output

clean_build:
	rm -rf build

clean_domains:
	rm -rf domains/*.pts.geo

clean_doxygen:
	rm -rf doxygen/html doxygen/latex doxygen/doxy_log.txt

clean_output:
	rm -rf output */output

clean_pycache:
	rm -rf __pycache__ */__pycache__ */*/__pycache__


## Generate the doxygen within the "doxygen" folder.
doxygen:
	cd doxygen; doxygen Doxyfile


## Format all .cxx and .hxx files that are parts of the library.
format:
	clang-format -i reproducables_python/parameters/*.hxx examples/parameters/*.hxx examples/*.cxx \
		tests_c++/*.cxx include/HyperHDG/*.hxx include/HyperHDG/*/*.hxx


## Download the submodules from GitHub.
submodules:
	git submodule update --init --recursive


## Perform a test that checks the library to work with all TEST_COMPILER.
test_all_compilers:
	$(foreach compiler, $(TEST_COMPILER), $(MAKE) CXX=${compiler} test_compiler;)

## Test the library to work with a certain compiler (given via the flag "comp").
test_compiler:
	$(MAKE) clean
	mkdir -p build	
	cd build; cmake -DCMAKE_CXX_COMPILER=$(CXX) ..
	cd build; make
	cd build; make test
