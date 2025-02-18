include config.mk

.PHONY: all-cpu util pscf fd1d pspc pspg test-cpu \
        clean clean-tests veryclean
# ======================================================================
# Main build targets

# Build all CPU and GPU programs
all:
	$(MAKE) all-cpu
	$(MAKE) pspg

# Build all programs that uses exclusively CPU hardware
all-cpu: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd fd1d; $(MAKE) all
	cd pspc; $(MAKE) all

# Build code in Util names (general scientific utilities)
util: 
	cd util; $(MAKE) all

# Build common code for PSCF in Util and Pscf namespaces
pscf: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all

# Build pscf_fd 1D finite difference program (install in BIN_DIR)
fd1d: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd fd1d; $(MAKE) all

# Build pscf_pcNd CPU code for periodic structures (install in BIN_DIR)
pspc: 
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd pspc; $(MAKE) all

# Build pscf_pcNd GPU code for periodic structures (install in BIN_DIR)
pspg:
	cd util; $(MAKE) all
	cd pscf; $(MAKE) all
	cd pspg; $(MAKE) all

# ======================================================================
# Test targets

# Run all CPU unit tests
test-cpu:
	cd util/tests; $(MAKE) all; $(MAKE) quiet
	cd pscf/tests; $(MAKE) all; $(MAKE) run
	cd fd1d/tests; $(MAKE) all; $(MAKE) run
	cd pspc/tests; $(MAKE) all; $(MAKE) run
	@cat util/tests/count > count
	@cat pscf/tests/count >> count
	@cat fd1d/tests/count >> count
	@cat pspc/tests/count >> count
	@echo " "
	@echo "Summary"
	@cat count
	@rm -f count

# ======================================================================
# Clean targets

# Remove object (*.o), dependency (*.d) and library (*.a) files
clean:
	cd util; $(MAKE) clean
	cd pscf; $(MAKE) clean
	cd fd1d; $(MAKE) clean

# Clean unit tests
clean-tests:
	cd util/tests; $(MAKE) clean
	cd pscf/tests; $(MAKE) clean
	cd fd1d/tests; $(MAKE) clean

# Remove all automatically generated files, recreate initial state
veryclean:
	cd util; $(MAKE) veryclean
	cd pscf; $(MAKE) veryclean
	cd fd1d; $(MAKE) veryclean
	cd pspc; $(MAKE) veryclean
	cd pspg; $(MAKE) veryclean
	rm -f util/config.mk
	rm -f pscf/config.mk
	rm -f fd1d/config.mk
	rm -f pspc/config.mk
	rm -f pspg/config.mk
ifneq ($(BLD_DIR),$(SRC_DIR))
	rm -f util/makefile
	rm -f pscf/makefile
	rm -f fd1d/makefile
	rm -f pspc/makefile
	rm -f pspg/makefile
	rm -f util/tests/makefile
	rm -f pscf/tests/makefile
	rm -f fd1d/tests/makefile
	rm -f pspc/tests/makefile
	rm -f configure
endif
	rm -f config.mk

# ======================================================================
