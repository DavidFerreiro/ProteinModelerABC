########################################################################
# - Compile ProteinEvolverABC -
########################################################################
# Makefile ProtABC
# 	Compile: make all
# 	Clean: make clean
#
PROGRAM_NAME = ProtModel


help:
	@echo ""
	@echo "To use the $(PROGRAM_NAME) Makefile, type"
	@echo "     make all           to compile the whole package but not install it"
	@echo "                          or remove the object files"
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo " "


all: 
	@echo "Creating bin folder ..."
	@mkdir bin
	@echo "Done!"
	@echo "Copying Gratham.csv file ..."
	@cp source/Gratham.csv ./bin
	@echo "Done!"
	@echo "Compiling DeltaGREM ..."
	$(MAKE) -C source/DeltaGREM_src clean
	$(MAKE) -C source/DeltaGREM_src all
	@cp source/DeltaGREM_src/DeltaGREM ./
	@echo "Done!"
	@echo ""
	@echo "Compiling ProteinEvolverProtABC ..."
	$(MAKE) -C source/ProteinEvolverProtABC clean
	$(MAKE) -C source/ProteinEvolverProtABC all
	@mkdir ./bin/Protein_EvolverABC
	@cp source/ProteinEvolverProtABC/ProteinEvolverProtABC1.2.0 bin/Protein_EvolverABC/
	@echo "Done!"
	@echo ""
	@echo "Note that R will require the following libraries: abc"
	@echo "These libraries can be installed from R by typing:"
	@echo "install.packages("abc")"
	@echo "And python will require the following libraries: os, sys, Biopython, rnamdom, numpy, warnings, pandas, csv, re, plataform and multiprocessing"
	@echo "These libraries can be installed from command line by typing:"
	@echo "pip install library_name"
	@echo "See the documentation for additional details about $(PROGRAM_NAME) and for cluster version information"
	@echo "Compilation completed!"
	@echo ""
	

clean:
	@echo "Removing executables ..."
	$(MAKE) -C source/DeltaGREM_src clean
	@rm -f DeltaGREM
	$(MAKE) -C source/ProteinEvolverProtABC clean
	@rm -r bin
	@echo "Done!"

