# differentiableImage
Computes the cubic interpolation reconstruction of an subsampled image function of the Gaussian blur variance of the original image

> Author: Guillaume Caron
>
> Dates: from February 2023 to ...

## Requirements : 

- [cmake](https://cmake.org/download/) (version 3.16.3 tested)
- [ceres-solver](https://github.com/ceres-solver/ceres-solver/tree/2.0.0) (version 2.0.0, github, master, commit 8d3e64dd5e64b346ed3e412cba2b70d760881bfb tested)
- [opencv](https://docs.opencv.org/4.x/d7/d9f/tutorial_linux_install.html) (version 4.2.0 tested)

## Usage

- To compile use:
     ```bash
     mkdir build 
     cd build
     cmake .. -DCMAKE_BUILD_TYPE=Release
     make 
     sudo make install 
     ```

- Use as library or run the program from the command line in the build directory:

```bash
./DifferentiableImageStandalone -desired {path_des} -current {path_cur}  {[optinals]-desiredDual {path_des_dual} -currentDual {path_cur_dual} -mask {path_mask} -display=0}
```
