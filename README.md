# flexWecDesignOpt
Automates WAMIT response analysis for use with wave energy converter structural design optimization

flexWecDesignOpt is a Python 3 package designed to automate hydrodynamic design analysis using the boundary element 
solver [WAMIT](http://wamit.com/). The package automates creating run configuration input files, writing .gdf file 
format mesh files, running WAMIT, and managing file structures to hold everything.

Usage requires creating a specifically formatted set of common input files and defining a custom class to define 
your design from a set of design parameters that can completely quantify your device. This device class is relatively 
quick to write for simple geometries, and is responsible for geometry creation and clarifying how to change the set of common input file values based off of the design parameters.

## Installation
flexWecDesignOpt can be installed with
```bash
pip install git+https://github.com/berriera/flexWecDesignOpt
```

## Features
1. Boundary element method input file text substitution
2. Mesh file creation and writing
3. Input file creation and file management
4. Additional tools for analyzing structurally flexible designs

## Dependencies
* Standard Python scientific libraries `numpy` and `scipy`
* `pygymsh` for parameterized surface meshing
* `meshmagick` for mesh trimming below the waterline and writing meshes to .gdf files

Note that `pygmsh` requires installing the application [Gmsh](http://gmsh.info/) and adding the application folder to your Python path


## Running in command line
All flexWecDesignOpt options can be run in the command line on built in device types using
```python
python main.py -i examples/WAMIT_meshing/input.yaml -m -r
```
##### Arguments
* `-i` or `--input` is a required argument and specifies the file path to the `input.yaml` file location
* Including `-m` or `--mesh` specifies to use the built-in meshing automation to write .gdf files based off the device geometry
* Including `-r` or `--run` specifies to run WAMIT on each set of input files to generate each set of output files

##### Creating a input.yaml file
The input .yaml file required for using the package in command line requires 4 values and has 3 optional ones.
1. Required values
    * `device_name`: the name of the device class to be used; current examples include a rigid-body `Cylinder` and
 structurally flexible `FlexibleBarge` class. You can use these example classes or add your own custom class to the file
 `device_types.py`
    * `common_file_directory`: the folder location of genericized WAMIT input files
    * `output_directory`: the folder location where all generated files will be stored
    * `cases_file`: the file name of the .csv file that lists the design variable quantities in the device class
2. Optional values:
    * `gmsh_exe_location`: location of the GMSH application folder used for surface mesh generation; should end 
    in `gmsh`
    * `mesh_refinement_factor`: a floating point number used for additional mesh refinement by GMSH; the default value 
    used in flexWecDesignOpt is 0.5 and values closer to 0.0 specify additional refinement
    * `run_wamit_command`: the command used to run WAMIT on your machine on a command line; on a Windows machine this 
    command is `C:\WAMTIv7\wamit`

## Use as a package
Refer to example file `flexible_barge_example.py` for an example of using flexWecDesignOpt in a script. The file uses 
Monte Carlo sampling to randomly test 10 different flexible barge designs and output them to the folder 
`examples/output`. The example can be run with the command
```bash
python flexible_barge_example.py
```
in the main package directory. Copying and pasting the code with the minor edits of removing commented out lines in the 
main loop and updating the top installation location variable for [Gmsh](http://gmsh.info/) will enable meshing and 
mesh file creation.

Note that example files have only been tested on Windows as of now.

### Creating a class for interfacing
Creating a wave energy converter design class requires making a substitution method and geometry method for translating
your design parameters into a full design. A simplified flexible barge example class is given below that uses the library 
[pygmsh](https://pypi.org/project/pygmsh/) for meshing. Refer to their 
[documentation](https://pygmsh.readthedocs.io/en/latest/index.html) for resources on how to create a parameterized geometry.

```python
import pygmsh
class FlexibleBarge(object):

    def __init__(self, design_vars):
        # These definitions translate how the length L, width w, and height h of the barge are translated from each row
        # of the .csv file
        self.L = design_vars[0]
        self.w = design_vars[1]
        self.h = design_vars[2]

        # Material constants can also be specified for the design or changed based on chosen values
        self.E = 30.720e6  # Modulus of elasticity, Pa
        self.rho = 500  # density, kg / (m ** 3)

        # Inertial relationships such as radii of gyration can be calculated using input values like this
        self.kx = ((1 / 12) * (self.L ** 2 + self.h ** 2)) ** (1 / 2)
        self.ky = ((1 / 12) * (self.w ** 2 + self.h ** 2)) ** (1 / 2)
        self.kz = ((1 / 12) * (self.L ** 2 + self.w ** 2)) ** (1 / 2)
        self.Cg = 0
        self.mass = self.rho * self.L * self.w * self.h

    def substitutions(self):
        # This method returns a dictionary that relates how text file substitutions correspond to device quantities
        return {'L': self.L, 'w': self.w, 'h': self.h, 
                'mass': self.mass, 'Cg': self.Cg, 
                'kx': self.kx, 'ky': self.ky, 'kz': self.kz
                }

    def geometry(self):
        # This method specifies how the device geometry relates to given input variables. For example, the barge is 
        # created as a box with dimensions L, w, and h centered at the origin
        geometry = pygmsh.opencascade.Geometry()
        geometry.add_box(x0=[-1 / 2 * self.L, -1 / 2 * self.w, -1 / 2 * self.h],
                         extents=[self.L, self.w, self.h])
        return geometry
```

### Creating generic input files for WAMIT
Input files should be made in the typical format for WAMIT in a common folder, with all constant variables specified as
is. All variables to be altered using the `device.substition()` method are labeled using the formatting 
`?dictionary_variable_name?` to specify the appropriate variable to plug in. A simple example for writing a .frc file
using this method and the above `FlexibleBarge` class example is
```text
BARGE DESIGN FILE: L = ?L?, w = ?w?, h = ?h?
1 2 2 1 0 0 0 0 0       IOPTN
?Cg? VCG
?kx? 	  0.000000 0.000000
0.000000 ?ky? 	   0.000000
0.000000 0.000000  ?kz?      XPRDCT
0           NBETAH
0           NFIELD
```