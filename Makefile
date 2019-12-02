.PHONY: doxygen run

default:
	mkdir -p output
	./setup.sh

doxygen:
	cd doxygen; doxygen Doxyfile

run:	
	python3 Executable.py
