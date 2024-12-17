# *ProteinModelerABC*
*ProteinModelerABC* is an evolutionary framework to estimate the best-fitting substitution model of protein evolution, among a set of complex substitution models that cannot be accommodated in likelihood functions, with approximate Bayesian computation (ABC). *ProteinModelerABC* can be run using the terminal or a graphical users interface in LinuxOS or MacOS siystems or in computer clusters based on LinuxOS.

## Reference

When using ProteinModelerABC, cite:

Ferreiro D, Branco C, Arenas M. 2024. Selection among site-dependent structurally constrained substitution models of protein evolution by approximate Bayesian computation. Bioinformatics, 40, 3, btae096. [https://doi.org/10.1093/bioinformatics/btae096](https://doi.org/10.1093/bioinformatics/btae096).

## Acknowledgments

This work was supported by the Spanish Ministry of Science and Innovation through the Grant [PID2019-107931GA-I00/AEI/10.13039/501100011033]. D.F. was funded by a fellowship from the Xunta de Galicia [ED481A-2020/192]. 

## Documentation and Support
For more information about using *ProteinModelerABC*, a manual is [attached](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Documentation).

If still in doubt, please do not hesitate to contact us (ferreirogarciadavid@gmail.com) for any question.

## Download
Click here to download the [computer](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC) or [cluster](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC_Cluster) version

## Install
For a first and rapid test of the framework, we recommend installing and running it with a rapid example (details below in the section "Fast examples ready to be run on local computer"), thus using an example placed in [Fast_examples](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Fast_examples).

To compile the program go to the main directory (which must contain all the files of the [computer](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC) or [cluster](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC_Cluster) versions) and type:
```
make install
```
*(Compilation tested using 9.3.0 GNU Fortran (GCC) version, the compiler should include the library "libg2c")*

To compile the GUI go to the main directory and type:
```
make install-GUI
```
Automatically, some Python libraries are installed but the user must install the abc R package and the Cluster modules if desired:

|Name	|Language	|Version |
|---------------|---------------|---------------|
|abc	|R	|All| |-|
|Biopython	|Python	|All|
|numpy	|Python	|All|
|pandas	|Python	|All|
|matplotlib	|Python	|All|
|mpi4py	|Python	|Cluster|
|PyInstaller	|Python	|GUI|


## Main commands to execute
An analysis can be executed from anywhere just placing there the main executable file after compilation (ProteinModelerABC.py) and all the corresponding input files (including structures.in).
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
python3 ProteinModelerABC_Cluster.py
```
And next execute "*launch_PMABC.sh*" file:
```
sbatch launch_PMABC.sh
```
## Fast examples ready to be run on local computer
The [first Fast-example](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Fast_examples/Fast_example1_Coalescent), [second Fast-example](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Fast_examples/Fast_example2_PhylogeneticTree) and [third_Fast-example](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Fast_examples/Fast_example3_ABCModels) are almost ready to be run. The user should go to the corresponding *Input* directory and:
```
make install
```
Next, the user must execute the framework typing:
```
python3 ProteinModelerABC.py
```
*(Note that the Settings.txt file is filled with no-biological meaning information)*


If the user prefers running the GUI must type:
```
make install
make install-GUI
```
Next, the user must execute the GUI:
```
Go to the "Executable" folder (created after GUI compilation) located inside GUI folder and clik the executable file
```

or

```
cd GUI/
python3 ProteinModelerABC_GUI.py
```
*(Note that the Settings.txt will be filled by the user)*

Each example includes an *Output* folder with an example of the output files. Note that to obtain the outputs under the three ABC methods the user must execute the ABCAnalysis.r file with the desire ABC method and with the appropiate ABC tolerance value (see [documentation](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Documentation)).

## Disclaimer
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version (at your option) of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place – Suite 330, Boston, MA 02111-1307, USA.
