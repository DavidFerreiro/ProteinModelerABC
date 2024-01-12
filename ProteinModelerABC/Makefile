########################################################################
# - Compile ProteinModelerABC -
########################################################################
#
PROGRAM_NAME = ProteinModelerABC
OS := $(shell uname)

CURDIR := $(shell pwd)

help:
	@echo ""
	@echo "To use the $(PROGRAM_NAME) Makefile, type"
	@echo "     make install          to compile the whole package "
	@echo "     make install-GUI          to compile the GUI "
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo " "


install: 
	@echo "Creating bin folder ..."
	@mkdir bin
	@echo "Done!"
	@echo "Copying Grantham.csv file ..."
	@cp source/Grantham.csv ./bin
	@echo "Done!"
	@echo "Copying Scripts folder ..."
	@cp -r source/Scripts ./bin
	@python bin/Scripts/Route.py $(CURDIR)
	@rm bin/Scripts/Route.py
	@echo "Done!"
	@echo "Compiling DeltaGREM ..."
	$(MAKE) -C source/DeltaGREM_src clean
	$(MAKE) -C source/DeltaGREM_src all
	@mkdir ./bin/DeltaGREM/
	@chmod +x source/DeltaGREM_src/DeltaGREM
	@cp source/DeltaGREM_src/DeltaGREM ./bin/DeltaGREM/
	@echo "Done!"
	@echo ""
	@echo "Compiling ProteinEvolverProtABC ..."
	$(MAKE) -C source/ProteinEvolverProtABC clean
	$(MAKE) -C source/ProteinEvolverProtABC all
	@mkdir ./bin/Protein_EvolverABC
	@chmod +x source/ProteinEvolverProtABC/ProteinEvolverProtABC1.2.0
	@cp source/ProteinEvolverProtABC/ProteinEvolverProtABC1.2.0 bin/Protein_EvolverABC/
	@chmod +x ProteinModelerABC.py
	@echo "Done!"
	@echo ""
	@echo ""
	@pip3 install Biopython
	@pip3 install numpy
	@pip3 install pandas
	@pip3 install matplotlib
	@echo "Succesfull intallation of Python libraries: Biopython, numpy, pandas and matplotlib"
	@echo "If any additional library is required the user should type:"
	@echo "pip install library_name"
	@echo "Note that R will require abc package"
	@echo "It can be installed from R by typing:"
	@echo "install.packages("abc")"
	@echo "See the documentation for additional details about $(PROGRAM_NAME) and for cluster version information"
	@echo "Compilation completed!"
	@echo ""

install-GUI:
	@pip3 install PyInstaller
	@echo "Creating Excutable folder ..."
	@mkdir ./GUI/Executable
	@echo "Done!"
	@echo "Creating ProteinModelerABC GUI executable ..."
	@python3 -m PyInstaller --onefile -w --add-data "./GUI/Mini.png:." --add-data "./GUI/Splash.png:."  ./GUI/ProteinModelerABC_GUI.py --exclude-module _bootlocale
ifeq ($(OS),Darwin)
	cp -r ./dist/ProteinModelerABC_GUI.app ./GUI/Executable 
else
	cp  ./dist/ProteinModelerABC_GUI ./GUI/Executable 
endif
	@rm -r dist
	@rm -r build
	@rm ProteinModelerABC_GUI.spec
	@echo "Compilation completed!"
	@echo ""
	@echo ""

clean:
	@echo "Removing executables ..."
	$(MAKE) -C source/DeltaGREM_src clean
	@rm -f DeltaGREM
	$(MAKE) -C source/ProteinEvolverProtABC clean
	@rm -r -f bin
	@rm -r -f ./GUI/Executable
	@echo "Done!"

