all: 
	$(MAKE) -C src all
	$(MAKE) -C STRIDE/src
	$(MAKE) -C TMalign/src


main: 
	$(MAKE) -C src


clean: 
	$(MAKE) -C src clean