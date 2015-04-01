# Example1
A example project that uses some parts of GLAMER

This simple example main() along with the sample parameter file will read in a fits density map and do various lensing calculations on it.

Once GLAMER has been installed it can be built in the same way using cmake.  The CMakeList.txt file is included.  
Linking you program to glamer is discussed in the wiki of the glenco/glamer project.

You will need a fits density map and you will need to change the name "map.fits" in the parameter file to this name.  The header requirements for this file are also given in the wiki of the glenco/glamer project.

For further information about the library functions see http://metcalf1.bo.astro.it
