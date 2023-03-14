# *ProteinModelerABC*
*ProteinModelerABC* is an evolutionary framework to estimate the best-fitting substitution model of protein evolution, among a set of complex substitution models that cannot be accommodated in likelihood functions, with approximate Bayesian computation (ABC). *ProteinModelerABC* can be run using the terminal or a graphical users interface in LinuxOS or MacOS siystems or in computer clusters based on LinuxOS.

## Documentation
For more information about using *ProteinModelerABC*, a manual is [attached](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Documentation).

## Download
Click here to download the [computer](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC) or [cluster](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC_Cluster) version

## Install
To compile the program go to the main directory and type:
```
make install
```
To compile the GUI go to the main directory and type:
```
make install-GUI
```
*(Note that the GUI compilation may fail with some Python versions)*

Automatically, some Python libraries are installed but the user must install the abc R package and the Cluster modules if desire:

|Name	|Language	|Version|
|---------------|---------------|---------------|
|abc	|R	|All|
|Biopython	|Python	|All|
|numpy	|Python	|All|
|pandas	|Python	|All|
|mpi4py	|Python	|Cluster|


## Main command to execute
1. To run it on command line you should type:
```
python3 ProteinModelerABC.py
```
2. To run it using GUI you should click on the executable file or type:
```
python3 ProteinModelerABC_GUI.py
```
3. To run it on cluster you should type:
```
python ProteinModelerABC_Cluster.py
```
And next execute "*launch_Simu.sh*" file:
```
sbatch launch_Simu.sh
```
## Fast examples ready to be run
The [Fast-examples](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Fast-examples) are almost ready to be run. The user should go to the directory and:
```
make install
```
Next, the user must execute the framework typing:
```
python3 ProteinModelerABC_GUI.py
```
*(Note that the Settings.txt file is filled with no-biological meaning information)*


If the user prefers running the GUI must type:
```
make install
make install-GUI
```
Next, the user must execute the GUI:

1.

Going to the "Executable" folder and cliking into the executable file 

or
   
2.
```
cd GUI/
python3 ProteinModelerABC_GUI.py
```
*(Note that the Settings.txt will be filled by the user)*

## Disclaimer
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version (at your option) of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place – Suite 330, Boston, MA 02111-1307, USA.
