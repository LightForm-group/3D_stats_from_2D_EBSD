# 3D_stats_from_2D_EBSD
MATALB script to extract 3D statistics from a 2D EBSD map read in with [MTEX](https://mtex-toolbox.github.io/download) to be read into [DREAM3D](http://dream3d.bluequartz.net/Download/).

## Equiaxed Processing
This can be used for EBSD maps where there is not a strong directionality in the grain axis. The program is currently set up to handle ctf data. There are a number of routines in the code that are summerized below. 

### Equivelent Diameter
Within DREAM3D the grain sizes are defined by a normal distribution in log(equivelent diameter). Here equivelent diameter is calculated via considering the volume of the grain as a sphere and finding its diameter. As a 2D slice intersects with the roughly sphereical grains this is a classical [Corpuscle Problem](https://pure.au.dk/portal/files/78832132/wicksell.pdf). This requires an adjustment ot the size distribution found in 2D using MteX using the [Saltykov algorithm](https://www.sciencedirect.com/science/article/pii/0026080073900013)

This leads to the following adjustment:
![image](https://user-images.githubusercontent.com/41897413/234549581-d894e686-9c25-4062-8c8e-33ee10a8fe22.png)

### Binning
![image](https://user-images.githubusercontent.com/41897413/234546528-641d5ca5-47cf-41ca-b9cf-121418c3248d.png)
*DREAM3D Size Distribution*

Within DREAM3D the statistics are split into bins based on grain size. This is done in the code by defining the same bin parameters and binning the adjusted grain equivelent diameters into the parameter 'diameters_binned' containing labels for each grain to the corresponding bin.

### Axis Ratios
For a 2D image this is challenging and currently it is assumed that the ratios are equal and that the distribution is equal in 2D and 3D. This is fitted to a Beta distribution. 

### Neighbor Distributions
This is achived by finding the number of neighbors in 2D within 1 grain diameter of each grain and then assuming a perfect 2D packing of circles to find distances between points which is generalised to perfect packing of spheres in 3D. The neighbors are then fitted to a lognormal distribution for each bin. 

### Axis ODF
This is challenging, the current method involves finding the distribution of angles of the long axis to the horizontal and assuming that the same distriubtution exists for all three euler angles giving a weight based on the relative probabilities. This cannot be fed into the GUI and must be fed into a DREAM3D script using an edited json file.

![image](https://user-images.githubusercontent.com/41897413/234552286-369e39ef-fbb2-4478-86ab-9692427be2a9.png)
*Distribution of angles in 2D*

### ODF and MDF
This comes easily from MTEX

## Fibrous Processing
This requires 2 EBSD maps which allows us to assume that the grains are roughly elongated boxes with two EBSD maps taken, one orthogonal to the long axis and one with the long axis in the plane. This gives us most of the structure but these maps must be related to each other. This is done by finding the distribution of the long axis from the NDRD map and then sampling from the distribution to assign a long axis value to the two short axis values found from the RDTD map.

This gives 3 axies and so 3D equivelent diameter, axis ratios and axis ODF are trivial. For the neighbor distributions we can assume that as the long axis will be longer than the equivelent diameter that all grain centres can be found in a 2D sampling of the RDTD map. 


