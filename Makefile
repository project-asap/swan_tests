TESTS=sparse_lu

SWANFLAGS=-fcilkplus -std=c++11
# SWANRTDIR= # supply on command line or in environment

all: test

sparse_lu: sparse_lu.cc
	$(CXX) $(SWANFLAGS) sparse_lu.cc -g -O2 -I$(SWANRTDIR)/include -L$(SWANRTDIR)/.libs -o sparse_lu

.PHONY: test
test: sparse_lu test.1 test.2 test.16

# MacOSX requires DYLD_LIBRARY_PATH to be set instead of LD_LIBRARY_PATH
#DYLD_LIBRARY_PATH=$(SWANRTDIR)/.libs:$(DYLD_LIBRARY_PATH) CILK_NWORKERS=$* ./sparse_lu > test.out.$* 2>&1

test.%: sparse_lu
	LD_LIBRARY_PATH=$(SWANRTDIR)/.libs:$(LD_LIBRARY_PATH) CILK_NWORKERS=$* ./sparse_lu > test.out.$* 2>&1
	@if grep "matrices are identical" test.out.$* > /dev/null ; then echo test with $* threads successful ; else echo test with $* threads failed ; fi


.PHONY: clean
clean:
	rm -fr sparse_lu sparse_lu.dSYM test.out.*
