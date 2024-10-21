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

- Create a new directory named build at the same level of differentiableImage directory
- Use cmake to fill the build directory :
     ```bash
     cd build
     cmake .. -DCMAKE_BUILD_TYPE=Release 
     ```
- Open the project or use the make command in the latter directory to build the exe file 
     ```bash
     make 
     sudo make install
     ```
- Run the programs from the command line from the differentiableImage parent directory, with arguments such as:

```bash
./DifferentiableImageStandalone -desired {path_des} -current {path_cur}  {[optinals] --display=false -mask mask}
```

- Plot residuals function of sigmas with gnuplot
```
shell: paste sigma.txt residuals.txt > sr.txt
```

gnuplot: 

plot 'sr.txt' using 1:2 with lines title "photometric reconstruction error", [0:30] 1 lc rgb "#FF0000" lw 2 title "1 diff intensity error", [0:30] 0.5 lc rgb "#00FF00" lw 2 title "0.5 diff intensity error"
set xrange [0:30]
set yrange [0:20]
replot

#set logscale y
