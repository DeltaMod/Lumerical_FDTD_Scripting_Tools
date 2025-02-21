# Lumerical_FDTD_Scripting_Tools
A collection of scripts categorised by their use-case intended for use in the Lumerical FDTD solver in the Ansys optics suite. Some scripts are cross-compatible with MODE as well.

Handling of exported .mat data is done using my other package called [DMU](https://github.com/DeltaMod/DMU/tree/main/DMU)
Each example setup file contains annoted python code that goes over how we use the package to analyse our data.

Generally speaking, to set up the environment to begin with we use:
```
from DMU import utils as dm
dm.Init_LDI() #To create data directories and settings json
#run dm.CUV(act="ddir") from console and select your data directories, then the next line loads the folders containing the data
filelist = dm.DataDir(act="load")["1"]
```
