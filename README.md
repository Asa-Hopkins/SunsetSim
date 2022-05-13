# SunsetSim
An attempt at simulating Rayleigh scattering. Place the sun anywhere in the sky to see what colours you get.

A fairly simple model is used where the Earth is spherical, the atmosphere is constant density and extends for 30km, and the sun's light contains equal amounts of red, green and blue.

Each frame is rendered by tracing a ray through each pixel, but rather than typical ray tracing, it instead calculates how much each particle it encounters is scattering light in that direction.

The result is therefore an integral, but calculating an integral for each pixel is prohibitative so a lot of approximations are used. The intention is to eventually make these approximations more reasonable, but even the current crude approximations look good.

To compare, numerical integration code is borrowed from https://www.codeproject.com/Articles/31550/Fast-Numerical-Integration, which is under the BSD license so I don't think there's any license issues.

## Controls
The controls are rather unintuitive, currently the mouse is used to look around and the arrow keys are used to move the sun. Left and right keys control the azimuthal angle, up and down keys control the polar angle. Press space to exit.

## Dependencies
Depends on CImg 2.9.2 or later.

## Installing

To install, clone the repository and set it to the current working directory. Then to compile for Linux run
```sh
g++ -O3 main.cpp -lpthread -lX11 -o Sunset
```
To compile for Windows, MinGW can be used, instead linking -lgdi32 instead of -lX11.

## Screenshots

Sunset            |  Daytime | Night
:-------------------------:|:-------------------------:|:-------------------------:
![](https://github.com/Asa-Hopkins/SunsetSim/blob/main/Screenshots/Screenshot1.png)  |  ![](https://github.com/Asa-Hopkins/SunsetSim/blob/main/Screenshots/Screenshot2.png) | ![](https://github.com/Asa-Hopkins/SunsetSim/blob/main/Screenshots/Screenshot3.png)
