PROJECT     	= HyperHDG
.PHONY:       	default clean distclean doxygen new run with_PythonCompileOptions object_files format \
								cython_cpp linking examples run_examples new_run_examples

make:
	$(MAKE) install
	
clean:
	$(MAKE) clean_build
	$(MAKE) clean_domains
	$(MAKE) clean_doxygen
	$(MAKE) clean_output
	$(MAKE) clean_pycache

clean_build:
	rm -rf build

clean_domains:
	rm -rf domains/*.pts.geo

clean_doxygen:
	cd doxygen; rm -rf html latex doxy_log.txt doxy_log.txt

clean_output:
	rm -rf output */output

clean_pycache:
	rm -rf __pycache__ */__pycache__

	rm -rf $(BUILD_DIR) $(OBJECT_DIR) $(CYTHON_DIR) $(CYTHON_FILE).c* $(DOXY_DIR) __pycache__
	rm -rf domains/*.pts.geo

doxygen:
	cd doxygen; doxygen Doxyfile

install:
	./setup.sh

submodules:
	git submodule update --init --recursive