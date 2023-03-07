########################################################################
# - Compile ProteinModelerABC GUI -
########################################################################
#
PROGRAM_NAME = ProteinModelerABC


help:
	@echo ""
	@echo "To create the $(PROGRAM_NAME) GUI executable file, type"
	@echo "     make all           to compile the whole package but not install it"
	@echo "                          or remove the object files"
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo " "


all: 
	@echo "Creating Excutable folder ..."
	@mkdir Executable
	@echo "Done!"
	@echo "Creating ProteinModelerABC GUI executable ..."
	@pip install pyinstaller
	@pyinstaller --onefile -w --add-data "Mini.png:." --add-data "Splash.png:."  ProteinModelerABC_GUI.py
	@cp -r ./dist/ProteinModelerABC_GUI.app ./Executable || cp -r ./dist/ProteinModelerABC_GUI ./Executable 
	@rm -r dist
	@rm -r build
	@rm ProteinModelerABC_GUI.spec
	@rm -r __pycache__
	@echo "Compilation completed!"
	@echo ""
	

clean:
	@echo "Removing executables ..."
	@rm -r Executable
	@echo "Done!"

