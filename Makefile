PROJECT     	= HDGonHYPERGRAPHS
.PHONY:       	default clean distclean doxygen new run with_PythonCompileOptions object_files format \
								cython_cpp linking examples run_examples new_run_examples

# Predefined directories for output, build files, and doxygen
SRC_DIR		= $(PWD)
OUTPUT_DIR	= $(SRC_DIR)/output
BUILD_DIR	= $(SRC_DIR)/build
DOXY_FILE_DIR	= $(SRC_DIR)/doxygen
EXAMPLE_DIR	= $(SRC_DIR)/tests_c++
INCLUDE_DIR = $(SRC_DIR)/include
PTCC_DIR = $(SRC_DIR)/submodules/tensor_product_chain_complex.git/include/

OBJECT_DIR  	= $(BUILD_DIR)/ObjectFiles
CYTHON_DIR  	= $(BUILD_DIR)/CythonFiles
EXAMPLE_BUILD	= $(BUILD_DIR)/tests_c++
CYTHON_FILE 	= ClassWrapper
DOXY_DIR	= $(DOXY_FILE_DIR)/html $(DOXY_FILE_DIR)/latex

# Extract relevant Python options, where overall version is chosen by user
PYTHON_VER		= 3
python_version_full  := $(wordlist 2,4,$(subst ., ,$(shell python$(PYTHON_VER) --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
python_version_patch := $(word 3,${python_version_full})

ifeq ($(PYTHON_VER), $(python_version_major))
	PYTHON      	= python$(python_version_major)
	PYTHON_M    	= /usr/include/python$(python_version_major).$(python_version_minor)
	CYTHONIZE			= cython
	CYTHONFLAGS		= -3 --cplus
else
$(error Python of version $(PYTHON_VER) needs to be installed on your computer!)
endif

# C++ Compiler options
COMPILER    	= g++
BASICFLAGS_R 	= -pthread -DNDEBUG -ggdb  -I$(PYTHON_M) -I$(SRC_DIR) -I$(INCLUDE_DIR) -I$(PTCC_DIR) \
								-fwrapv -O2 -Wall -g -fstack-protector-strong -Wformat -Werror=format-security \
								-Wdate-time -D_FORTIFY_SOURCE=2 -fPIC --std=c++17
BASICFLAGS  	= -pthread -ggdb  -I$(PYTHON_M) -I$(SRC_DIR) -I$(INCLUDE_DIR) -I$(PTCC_DIR) \
								-fwrapv -O2 -Wall -g -fstack-protector-strong -Wformat -Werror=format-security \
								-Wdate-time -D_FORTIFY_SOURCE=2 -fPIC --std=c++17
# Some additional useful flags: -Wpedantic -Wextra \

# C++ Linker options
LINKER      		= g++
LINKERPREFLAGS  = -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions \
                  -Wl,-z,relro -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong \
                  -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2
LINKERPOSTFLAGS = -llapack -lstdc++fs

# Sets of source and object files
SOURCE_FILES  := $(foreach src_dir, $(SRC_DIR), $(wildcard *.cxx))
OBJECTS       := $(foreach src, $(SOURCE_FILES), $(OBJECT_DIR)/$(src:.cxx=.o))
EXAMPLE_FILES	:= $(foreach src, $(EXAMPLE_DIR), $(wildcard $(EXAMPLE_DIR)/*.cxx))
EXAMPLE_HELP	:= $(foreach src, $(EXAMPLE_FILES), $(src:.cxx=.e))
EXAMPLE_OBJS	:= $(foreach src, $(EXAMPLE_HELP), $(subst $(EXAMPLE_DIR),$(EXAMPLE_BUILD),$(src)))
EXAMPLE_EXES	:= $(foreach src, $(EXAMPLE_OBJS), $(src:.e=.exe))
TEST_EXES			:= $(foreach src, $(EXAMPLE_BUILD), $(wildcard $(EXAMPLE_BUILD)/*.exe))

make:
	$(MAKE) run
	
clean:
	rm -rf $(BUILD_DIR) $(OBJECT_DIR) $(CYTHON_DIR) $(CYTHON_FILE).c* $(DOXY_DIR) __pycache__
	rm -rf domains/*.pts.geo

distclean:
	$(MAKE) clean
	rm -rf $(OUTPUT_DIR)

doxygen:
	cd $(DOXY_FILE_DIR); doxygen Doxyfile

new:
	$(MAKE) clean
	make

run:
	$(PYTHON) executable.py

format:
	clang-format -i tests_c++/*.cxx cython/*.cxx include/HyperHDG/*.hxx include/HyperHDG/*/*.hxx

tests:
	mkdir -p $(EXAMPLE_BUILD)
	$(MAKE) example_objects
	$(MAKE) example_linking

run_tests:
	$(MAKE) tests
	./build/tests_c++/plot_simple.exe;
	./build/tests_c++/plot_refined.exe;
	./build/tests_c++/point_functions.exe;
	./build/tests_c++/diffusion_floatings.exe;
	./build/tests_c++/diffusion_plotting.exe;
	./build/tests_c++/diffusion_dimensions.exe;
	./build/tests_c++/elasticity_simple.exe;
	$(PYTHON) tests_python/diffusion_uniform.py

new_run_tests:
	$(MAKE) clean
	$(MAKE) run_tests


example_objects: $(EXAMPLE_OBJS)

$(EXAMPLE_BUILD)/%.e: $(EXAMPLE_DIR)/%.cxx
	$(COMPILER) $(BASICFLAGS) -c $^ -o $@

example_linking: $(EXAMPLE_EXES)

$(EXAMPLE_BUILD)/%.exe: $(EXAMPLE_BUILD)/%.e
	$(LINKER) $^ -o $@ $(LINKERPOSTFLAGS)


print_variables:
	@echo 'SRC_DIR=    ' $(SRC_DIR)
	@echo 'OBJECT_DIR= ' $(OBJECT_DIR)
