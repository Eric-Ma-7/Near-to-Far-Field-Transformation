# Dirac Vortex Topological Cavities
Research on the optical far field of Dirac vortex topological cavities

## Near to far field transformation

### Introduction

MATLAB class `Field2D` supports the some operations on near to far field transformation.

There are a few public properties in `Field2D`:
+ Properties `x`, `y`, `nx`, `ny` are matrices of the same size to store the data of near field. `x` and `y` are the coordinate mesh while `nx` and `ny` are the x- and y-components of the near field.
+ Properties `kx`, `ky`, `fx`, `fy` are matrices of the same size to store the data of far field. `fx` and `fy` are x- and y-components of far field .Different with near field, the coordinate of far field is indicated by the wave vector `kx` and `ky`.
+ Properties `freq0`, `lambda0`, `k0` are the frequency, wavelength, wave vectore in vacuum respectively.

### Methods Declaration

```matlab
obj = Field2D(x_in, y_in, fx_in, fy_in, freq, unit)
```
The constructor accepts the input data of near field to assign a class.
+ `x_in`: the x-coordinate of near field. Must be a n-vector.
+ `y_in`: the x-coordinate of near field. Must be a n-vector.
+ `fx_in`: the x-component of near field. Must be a n-vector.
+ `fy_in`: the x-component of near field. Must be a n-vector.
+ `freq`: the frequency of the input near field. Must be a real number.
+ `unit`: the unit of x- and y-coordinate. Valid values are `'m', 'cm', 'mm', 'um', 'nm', 'a'`.

```matlab
obj = obj.getFarFieldCPU({kx, ky})
obj = obj.getFarFieldCPU({kx, ky}, CoreNumber)
obj = obj.getFarFieldCPU(farGridNumber)
obj = obj.getFarFieldCPU(farGridNumber, CoreNumber)
```
The method `getFarFieldCPU(-)` uses CPU to calculate parallelly and save the far field of the near field store in class `obj`.
+ `{kx, ky}`: this is a single input variable, which is a cell in MATLAB. The first element is a n-vector contains the value of `kx` and the second element is the same.
+ `farGridNumber`: the number of grids of far field. It is no different with `kx=ky=linspace(-1,1,farGridNumber)`.
+ `CoreNumber`: the CPU core number of parallel computation. Default is the half of total CPU core.

This method will first check if a parallel pool is running or not. if not, this method will start a parallel pool with size of `CoreNumber`.

```matlab
obj = obj.getFarFieldGPU(fargrid)
obj = obj.getFarFieldGPU({kx, ky})
```
The method `getFarFieldGPU(-)` uses GPU to compute and will return the same result as `getFarFieldCPU(-)`. `getFarFieldGPU(-)` is recommended to use because of its high efficiency.

```matlab
obj.plotNearField()
obj.plotNearField(colormap)
s = obj.plotNearField(-)
```
The method `plotNearField(-)` plots the near field power distribution `abs(obj.nx).^2+abs(obj.ny).^2`, and returns a image object of class `Patch`.
+ `colormap`: the colormap. Default is `parula`

```matlab
obj.plotFarField(quiverSize,density,range)
s = obj.plotFarField(-)
```
The method `plotFarField(-)` plots the far field power distribution `P=abs(obj.fx).^2+abs(obj.fy).^2`, and  returns a image object of class `Patch`. This method also plots the polarization of the far field in green line.
+ `quiverSize`: this factor determines size of polarization marks (green line). `quiverSize=1` is recommended.
+ `density`: this factor determines the density of polarization marks. This factor must be an interger (1,2,3...). A greater `density` indicates a smaller marker density. `density=2` is recommended
+ `range`: this factor determines the ploting range of polarization marks. The marks will be plotted in the area where the internsity `P<Pmax*range`. A smaller `range` indicates a bigger marked area. `range=0.1` is recommended.

### How to use

First, load data from file and use the constructor to assign an object `field`
```matlab
clear
data = readmatrix('near-field-data.csv');
x_in = data(:,1);
y_in = data(:,1);
Ex_in = real(data(:,3));
Ey_in = real(data(:,4));
f0 = 3.0454e14;
field = Field2D(x_in, y_in, Ex_in, Ey_in, f0, 'nm');
```
Then use `getFarFieldGPU(-)` to calculate the far field
```matlab
field = field.getFarFieldGPU(101);
```
Finally, use `plotFarField(-)` to plot the far field power distribution.
```matlab
figure
s = field.plotFarField(1,2,0.1);
```

***

This project is licensed under [GPL v3.0 license](https://www.gnu.org/licenses/gpl-3.0.en.html)
