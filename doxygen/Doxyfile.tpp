# The documentation of the submodule tensor_product_polynomials, linked from HyperHDG's
# through its tag file. Run from the source root before doxygen/Doxyfile, with
# DOXYGEN_OUTPUT_DIRECTORY set (see there); the pages go to $(DOXYGEN_OUTPUT_DIRECTORY)/html/tpp.

PROJECT_NAME         = "TPP — tensor product polynomials"
OUTPUT_DIRECTORY     = $(DOXYGEN_OUTPUT_DIRECTORY)
HTML_OUTPUT          = html/tpp
GENERATE_TAGFILE     = $(DOXYGEN_OUTPUT_DIRECTORY)/tpp.tag
INPUT                = submodules/tensor_product_polynomials.git/include
RECURSIVE            = YES

TAB_SIZE             = 2
EXTRACT_ALL          = YES
EXTRACT_PRIVATE      = YES
EXTRACT_STATIC       = YES
USE_MATHJAX          = YES
GENERATE_LATEX       = NO
HAVE_DOT             = YES
CALL_GRAPH           = YES
