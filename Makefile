PROJECT     	= HDGonHYPERGRAPHS
.PHONY:       	default clean distclean doxygen new run with_PythonCompileOptions object_files \
								cython_cpp linking

SRC_DIR     	= .
OUTPUT_DIR		= output
BUILD_DIR   	= build
DOXY_FILE_DIR	= doxygen

OBJECT_DIR  	= $(BUILD_DIR)/ObjectFiles
CYTHON_DIR  	= $(BUILD_DIR)/CythonFiles
CYTHON_FILE 	= ClassWrapper
DOXY_DIR			= $(DOXY_FILE_DIR)/html $(DOXY_FILE_DIR)/latex

COMPILER    	= g++
BASICFLAGS  	= -pthread -DNDEBUG -g -fwrapv -O2 -Wall -g -fstack-protector-strong -Wformat \
								-Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -fPIC \
								--std=c++17 -I$(PYTHON_M)

LINKER      		= g++
LINKERPREFLAGS  = -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions -Wl,-Bsymbolic-functions \
                  -Wl,-z,relro -Wl,-Bsymbolic-functions -Wl,-z,relro -g -fstack-protector-strong \
                  -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2
LINKERPOSTFLAGS = -llapack

PYTHON      	= python3m
PYTHON_M    	= /usr/include/python3.6m
CYTHONIZE			= cython
CYTHONFLAGS		= -3 --cplus


SOURCE_FILES  := $(foreach src_dir, $(SRC_DIR), $(wildcard *.C))
OBJECTS       := $(foreach src, $(SOURCE_FILES), $(OBJECT_DIR)/$(src:.C=.o))


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


object_files: $(OBJECTS)

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.C
	$(COMPILER) $(BASICFLAGS) -c $^ -o $@

cython_cpp: $(CYTHON_DIR)/$(CYTHON_FILE).o

$(CYTHON_DIR)/%.o: $(SRC_DIR)/$(CYTHON_FILE).pyx $(SRC_DIR)/$(CYTHON_FILE).pxd2
	rm -rf $(CYTHON_DIR)/*
	cp $(SRC_DIR)/$(CYTHON_FILE).pyx $(CYTHON_DIR)/$(CYTHON_FILE).pyx
	cp $(SRC_DIR)/$(CYTHON_FILE).pxd2 $(CYTHON_DIR)/$(CYTHON_FILE).pxd
	cd $(CYTHON_DIR); $(CYTHONIZE) $(CYTHONFLAGS) $(CYTHON_FILE).pyx
	cd $(CYTHON_DIR); $(COMPILER) $(BASICFLAGS) -c $(CYTHON_FILE).cpp -o $(CYTHON_FILE).o

linking: $(BUILD_DIR)/$(CYTHON_FILE).so

$(BUILD_DIR)/%.so: $(OBJECT_DIR)/*.o $(CYTHON_DIR)/*.o
	$(LINKER) $(LINKERPREFLAGS) $^ -o $@ $(LINKERPOSTFLAGS)
