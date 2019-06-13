# Painting Analysis

This is an interface to find vanishing points and lines in the image of a painting. The interface use [Qt](https://www.qt.io/). Then you can simply compile and run the code with the following commands:
```bash
mkdir build
cd build
qmake ..
make
./interface
```

## Requirements

The code is Parallelized with [OpenMP](https://www.openmp.org/). On a Debian system you can run the following command to install it:
```bash
sudo apt install libomp-dev
```