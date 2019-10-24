BINDIR=$(abspath "./bin")
export BINDIR

all: 
	$(MAKE) -C src/ProteinEvaluation all
	
	if [ ! -d "src/STRIDE/src" ]; then tar -xzf "src/STRIDE/stride.tar" -C "src/STRIDE"; fi
	$(MAKE) -C src/STRIDE/src
	
	$(MAKE) -C src/TMalign


main: 
	$(MAKE) -C src/ProteinEvaluation all


clean: 
	$(MAKE) -C src/ProteinEvaluation clean