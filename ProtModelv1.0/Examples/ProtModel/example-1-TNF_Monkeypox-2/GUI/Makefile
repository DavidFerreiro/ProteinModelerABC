########################################################################
# - Compile ProtModel -
########################################################################
#
PROGRAM_NAME = ProtModel


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
	@echo "Creating ProtModel GUI executable ..."
	@pip install pyinstaller
	@pyinstaller --onefile -w --add-data "Mini.png:." --add-data "Splash.png:."  ProtModel_GUI.py
	@cp -r ./dist/ProtModel_GUI.app ./Executable
	@rm -r dist
	@rm -r build
	@rm ProtModel_GUI.spec
	@echo "Compilation completed!"
	@echo ""
	

clean:
	@echo "Removing executables ..."
	$(MAKE) -C source/DeltaGREM_src clean
	@rm -r Executables
	@echo "Done!"

