all: 
	$(MAKE) -C src
	$(MAKE) -C STRIDE/src
	$(MAKE) -C TMalign/src


main: 
	$(MAKE) -C src


clean: 
	$(MAKE) -C src clean