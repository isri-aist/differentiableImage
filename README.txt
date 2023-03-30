##########################################################################################
#
# differentiableImage: computes the cubic interpolation reconstruction of an subsampled image function of the Gaussian blur variance of the original image
#
# Author: Guillaume Caron
#
# Dates: from February 2023 to ...
##########################################################################################

0. Prerequisites to run these codes: 
     install cmake (version 3.16.3 tested)
     install ceres-solver (version 2.0.0, github, master, commit 8d3e64dd5e64b346ed3e412cba2b70d760881bfb tested)
     install opencv (version 4.2.0 tested)

1. create a new directory named differentiableImage-build at the same level of differentiableImage directory
2. use cmake to fill the build directory in: cd differentiableImage-build; cmake ../differentiableImage/ -DCMAKE_BUILD_TYPE=Release 
3. open the project or use the make command in the latter directory to build the exe file
4. run the programs from the command line from the differentiableImage parent directory, with arguments such as:

	./differentiableImage	-desired path/to/image.jpg

Options:
-isEqui true: if the image considered is equirectangular

