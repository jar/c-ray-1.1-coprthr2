# C-Ray 1.1 for Parallella!
By James A Ross

![C-Ray rendering on Parallella Epiphany-III cores](foo.png)
 
## Original Code
http://www.futuretech.blinkenlights.nl/c-ray.html

## Prerequisites

* [COPRTHR-2 SDK](http://www.browndeertechnology.com/coprthr2.htm) for host and coprocessor code
* [SDL 2.0](https://www.libsdl.org/download-2.0.php) (optional) for SDL real-time rendering view

## Usage
```
    compile:  just type 'make'. 
...with SDL:  'USESDL=1 make'
        run:  ./c-ray-mt
    options:   -t <num>   how many threads to use (default: 1)
               -s WxH     where W is the width and H the height of the image
               -r <rays>  shoot <rays> rays per pixel (antialiasing)
               -i <file>  read from <file> instead of stdin
               -o <file>  write to <file> instead of stdout
               -h         help
```
Alternatively, the following syntax for scene input and image
output are equivalent:
```
cat scene | ./c-ray-mt -t 16 > foo.ppm
./c-ray-mt -t 16 -o foo.ppm -i scene
```

## Input Scene Format
```
Scene file format:
  # sphere (many)
  s  x y z  rad   r g b   shininess   reflectivity
  # light (many)
  l  x y z
  # camera (one)
  c  x y z  fov   tx ty tz
```

## Benchmarks
| Command                                                            | Execution Time (ms) |
|:------------------------------------------------------------------:|:-------------------:|
| <code>./c-ray-mt -t 16 -s 800x600 -r 1  -o foo.ppm -i scene</code> |                 318 |
| <code>./c-ray-mt -t 1  -s 800x600 -r 1  -o foo.ppm -i scene</code> |                5069 |
| <code>./c-ray-mt -t 16 -s 800x600 -r 16 -o foo.ppm -i scene</code> |                4988 |

By comparison, the [original C-Ray
1.1](http://www.futuretech.blinkenlights.nl/c-ray.html) code executed the first
benchmark on the dual ARM Cortex-A9 cores (with four threads) in 1360 ms, 4.3x
slower.

## Major Modifications from Original C-Ray 1.1
As Epiphany-III is a coprocessor some modifications had to be made:
* [Single Precision](#single_precision)
* [Host Code for Epiphany Offload](#host_code)
* [Parallel Load Balancing](#load_balancing)
* [Miscellaneous Supporting Code](#supporting_code)
* [Image Output](#image_output)
 
### <a name="single_precision"></a>Single Precision (32-bit)
The original code used 64-bit double precision instructions in the ray tracer.
The Epiphany-III cores within the Parallella board do not have hardware support
for double precision so it's not a completely fair comparison.

### <a name="host_code"></a>Host Code for Epiphany Offload
The host and device code use the [COPRTHR-2
SDK](http://www.browndeertechnology.com/coprthr2.htm).  Unified Virtual Memory
(UVA) is used to simplify sharing of data between the ARM host and Epiphany
coprocessor.  There are no explicit memory copies on the host with the
exception of marshalling data.  This is done to simplify the copy from shared
DRAM to local SRAM within each Epiphany core. 

### <a name="load_balancing"></a>Parallel Load Balancing
Parallel work is typically distributed evenly across multiple cores, but this
is not ideal for applications which have data-dependent computation -- some
cores would finish faster than others and go idle. One solution to this problem
is to use a mutex to fetch and increment a counter located within one of the
cores.  The relatively minor overhead of this operation keeps all cores busy
for the duration of the computation.

### <a name="supporting_code"></a>Miscellaneous Supporting Code
Several supporting math routines are used rather than the defaults from the
compiler.  The routines, `_rsqrt`, `_sqrt`, and `_inv` calculate the inverse
square root, square root, and inverse, respectively. They presently use four
newton iterations for higher precision. There are three other supporting math
routines written by SUN: `__powf`, `__copysignf`, and `__scalbnf` that are
included within the code. These routines are likely targets for optimization.

### <a name="image_output"></a>Image Output
A real-time viewer has been added with the SDL2 interface to visualizing the
state of the rendering. It presently uses a 5Hz update rate and is most useful
for long-running jobs. You can see the effect of the parallel load balancing. A
separate thread is created for the display loop.  In order to enable SDL, you
must have a clean build and define the USESDL value.
```
make distclean
USESDL=1 make
# An example of a long-running visualization:
./c-ray-mt -t 16 -r 128 -o foo.ppm -i scene
```

## Known Issues
Likely due to the conversion to 32-bit single precision floating point, some
residual graphical errors (so called 'surface acne') may appear. One can adjust
the error margin value, `ERR_MARGIN`, within the code, which may improve the
result.
