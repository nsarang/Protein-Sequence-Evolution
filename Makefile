all: 
	$(MAKE) -C src all
	
	if [ ! -d "./STRIDE/src" ]; then tar -xzf "./STRIDE/stride.tar" -C "./STRIDE"; fi
	$(MAKE) -C STRIDE/src
	
	$(MAKE) -C TMalign/src


main: 
	$(MAKE) -C src


clean: 
	$(MAKE) -C src clean