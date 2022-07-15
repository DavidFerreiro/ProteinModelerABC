# *ProtModel*
*ProtModel* is a framework for selecting the best-fitting substitution model of evolution for protein alignments considering both empirical and structuraly constrained substitution models by approximate Bayesian computation. *ProtModel* can be run using the terminal or a graphical users interface in Linux or MacOS siystems or in clusters.

## Documentation
For more information about using *ProtModel*, a manual is [attached](https://github.com/DavidFerreiro/ProtModel/tree/main/ProtModelv1.0/Documentation).

## Citation
If you use *ProtModel*, please cite the following:

## Download
Click here download the [computer](https://github.com/DavidFerreiro/ProtModel/tree/main/ProtModelv1.0/ProtModel) or [cluster](https://github.com/DavidFerreiro/ProtModel/tree/main/ProtModelv1.0/ProtModel_Cluster) version

## Install
To install the DeltaGREM and ProteinEvolver frameworks you should simply run:
```
make all
```
## Execute
1. To run it on command line you should type:
```
python3.9 ProtModelGeneral-M.py
```
2. To run it using GUI you should click on the executable file or type::
```
python3.9 ProtModel_GUI.py
```
3. To run it on cluster you should type:
```
python ProtModel_Cluster.py
```
And next execute "*launch_simu.sh*" file:
```
sbatch launch_simu.sh
```
## Disclaimer
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version (at your option) of the License. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place – Suite 330, Boston, MA 02111-1307, USA.
