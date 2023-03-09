# *ProteinModelerABC*
*ProteinModelerABC* is an evolutionary framework to estimate the best-fitting substitution model of protein evolution, among a set of complex substitution models that cannot be accommodated in likelihood functions, with approximate Bayesian computation (ABC). *ProteinModelerABC* can be run using the terminal or a graphical users interface in LinuxOS or MacOS siystems or in computer clusters based on LinuxOS.

## Documentation
For more information about using *ProteinModelerABC*, a manual is [attached](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/Documentation).

## Download
Click here to download the [computer](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC) or [cluster](https://github.com/DavidFerreiro/ProteinModelerABC/tree/main/ProteinModelerABC_Cluster) version

## Install
To install the DeltaGREM and ProteinEvolver frameworks and compile the graphical user interface (GUI) you should simply run:
```
make all
```
Next, Python3, R and some libraries and modules have to be installed:

|Name	|Language	|Version|
|---------------|---------------|---------------|
|abc	R	|All|
|os	Python	|All|
|sys	Python	|All|
|Biopython	Python	1All|
|random	Python	|All|
|numpy	Python	|All|
|warnings	Python	|All|
|pandas	Python	|All|
|csv	Python	|All|
|multiprocessing	Python	|Command line and GUI|
|re	Python	|All|
|platform	Python	|All|
|mpi4py	Python	|Cluster|
|threading	Python	|GUI|
|tkinter	Python	|GUI|
|time	|Python	|GUI|

## Execute
1. To run it on command line you should type:
```
python3.9 ProteinModelerABC.py
```
2. To run it using GUI you should click on the executable file or type:
```
python3.9 ProteinModelerABC_GUI.py
```
3. To run it on cluster you should type:
```
python ProteinModelerABC_Cluster.py
```
And next execute "*launch_Simu.sh*" file:
```
sbatch launch_Simu.sh
```
## Disclaimer
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version (at your option) of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place – Suite 330, Boston, MA 02111-1307, USA.
