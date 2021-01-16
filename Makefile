PROJECT     	= HyperHDG
.PHONY:       	build clean distclean clean_build clean_domains clean_doxygen clean_output \
								clean_pycache doxygen format submodules test_all_compilers test_github test_compiler


TEST_COMPILER = clang++-9 clang++-10 g++-9 g++-10


make:
	$(MAKE) build


build:
	mkdir -p build
	cd build; cmake ..
	cd build; make
	cd build; make test


clean:
	$(MAKE) clean_build
	$(MAKE) clean_domains
	$(MAKE) clean_doxygen
	$(MAKE) clean_pycache

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


doxygen:
	cd doxygen; doxygen Doxyfile


format:
	clang-format -i reproducables_python/parameters/*.hxx examples/parameters/*.hxx examples/*.cxx \
		tests_c++/*.cxx include/HyperHDG/*.hxx include/HyperHDG/*/*.hxx


submodules:
	git submodule update --init --recursive


test_all_compilers:
	$(foreach compiler, $(TEST_COMPILER), $(MAKE) test_compiler comp=${compiler};)

test_github:
	$(MAKE) clean
	mkdir -p build
	cd build; cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_CXX_FLAGS="-DNOFILEOUT" \
		-DNOPYTHONTESTS=True ..
	cd build; make
	cd build; make test

test_compiler:
	$(MAKE) clean
	mkdir -p build	
	cd build; cmake -DCMAKE_CXX_COMPILER=$(comp) ..
	cd build; make
	cd build; make test