PROJECT     	= HDGonHYPERGRAPHS
.PHONY:       	default clean distclean doxygen new run with_PythonCompileOptions object_files \
								cython_cpp linking examples run_examples new_run_examples

# Predefined directories for output, build files, and doxygen
SRC_DIR     	= .
OUTPUT_DIR		= output
BUILD_DIR   	= build
DOXY_FILE_DIR	= doxygen
EXAMPLE_DIR		= examples_c++

OBJECT_DIR  	= $(BUILD_DIR)/ObjectFiles
CYTHON_DIR  	= $(BUILD_DIR)/CythonFiles
EXAMPLE_BUILD	= $(BUILD_DIR)/C++ExampleBuild
CYTHON_FILE 	= ClassWrapper
DOXY_DIR			= $(DOXY_FILE_DIR)/html $(DOXY_FILE_DIR)/latex

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
BASICFLAGS  	= -pthread -DNDEBUG -g -fwrapv -O2 -Wall -g -fstack-protector-strong -Wformat \
								-Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -fPIC \
								--std=c++17 -I$(PYTHON_M)

# C++ Linker options
LINKER      		= g++
LINKERPREFLAGS  = -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions \
                  -Wl,-z,relro -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong \
                  -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2
LINKERPOSTFLAGS = -llapack

# Sets of source and object files
SOURCE_FILES  := $(foreach src_dir, $(SRC_DIR), $(wildcard *.C))
OBJECTS       := $(foreach src, $(SOURCE_FILES), $(OBJECT_DIR)/$(src:.C=.o))
EXAMPLE_FILES	:= $(foreach src, $(EXAMPLE_DIR), $(wildcard $(EXAMPLE_DIR)/*.C))
EXAMPLE_HELP	:= $(foreach src, $(EXAMPLE_FILES), $(src:.C=.e))
EXAMPLE_OBJS	:= $(foreach src, $(EXAMPLE_HELP), $(subst $(EXAMPLE_DIR),$(EXAMPLE_BUILD),$(src)))
EXAMPLE_EXES	:= $(foreach src, $(EXAMPLE_OBJS), $(src:.e=.exe))
TEST_EXES			:= $(foreach src, $(EXAMPLE_BUILD), $(wildcard $(EXAMPLE_BUILD)/*.exe))


default:
	mkdir -p $(OUTPUT_DIR) $(BUILD_DIR) $(OBJECT_DIR) $(CYTHON_DIR)
	make object_files
	make cython_cpp
	make linking

clean:
	rm -rf $(BUILD_DIR) $(OBJECT_DIR) $(CYTHON_DIR) $(CYTHON_FILE).c* $(DOXY_DIR)

distclean:
	make clean
	rm -rf $(OUTPUT_DIR)

doxygen:
	cd $(DOXY_FILE_DIR); doxygen Doxyfile

new:
	make clean
	make

run:
	make
	PYTHONPATH=$(BUILD_DIR) $(PYTHON) Executable.py

with_PythonCompileOptions:
	mkdir -p $(OUTPUT_DIR)
	$(PYTHON) PythonCompileOptions.py build_ext --inplace

examples:
	make
	mkdir -p $(EXAMPLE_BUILD)
	make example_objects
	make example_linking

run_examples:
	make examples
	./build/C++ExampleBuild/DiffusionTest1.exe;
	./build/C++ExampleBuild/DiffusionTest2.exe;

new_run_examples:
	make clean
	make run_examples


example_objects: $(EXAMPLE_OBJS)

$(EXAMPLE_BUILD)/%.e: $(EXAMPLE_DIR)/%.C
	$(COMPILER) --std=c++17 -c $^ -o $@

example_linking: $(EXAMPLE_EXES)

$(EXAMPLE_BUILD)/%.exe: $(OBJECT_DIR)/*.o $(EXAMPLE_BUILD)/%.e
	$(LINKER) $^ -o $@ $(LINKERPOSTFLAGS)

object_files: $(OBJECTS)

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.C
	$(COMPILER) $(BASICFLAGS) -c $^ -o $@

cython_cpp: $(CYTHON_DIR)/$(CYTHON_FILE).o

$(CYTHON_DIR)/%.o: $(SRC_DIR)/$(CYTHON_FILE).pyx $(SRC_DIR)/$(CYTHON_FILE).pxd2 $(SRC_DIR)/*.hpp
	rm -rf $(CYTHON_DIR)/*
	cp $(SRC_DIR)/$(CYTHON_FILE).pyx $(CYTHON_DIR)/$(CYTHON_FILE).pyx
	cp $(SRC_DIR)/$(CYTHON_FILE).pxd2 $(CYTHON_DIR)/$(CYTHON_FILE).pxd
	cd $(CYTHON_DIR); $(CYTHONIZE) $(CYTHONFLAGS) $(CYTHON_FILE).pyx
	cd $(CYTHON_DIR); $(COMPILER) $(BASICFLAGS) -c $(CYTHON_FILE).cpp -o $(CYTHON_FILE).o

linking: $(BUILD_DIR)/$(CYTHON_FILE).so

$(BUILD_DIR)/%.so: $(OBJECT_DIR)/*.o $(CYTHON_DIR)/*.o
	$(LINKER) $(LINKERPREFLAGS) $^ -o $@ $(LINKERPOSTFLAGS)
