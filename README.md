### DICCOR ###

A full set of MATLAB routines for generating displacement vector fields from time series of moving textures. The algorithm uses a spatial Fourier transform to efficiently compute the cross-correlation function for two images, as well as a variety of spatial and temporal smothing methods to clean up the resulting time-resolved vector field.

The code was developed by the Brangwynne laboratory at Princeton University. If using this code please cite:

     W.Gilpin, S. Uppaluri, C. Brangwynne “Worms under pressure: 
     bulk mechanical properties of C. elegans are independent of the cuticle” 
     Biophysical Journal, 2015.

Sample Displacement Field        |  Sample Maximal Strain Projection
:-------------------------:|:-------------------------:
![](sample_output/vec_field/overlay4.png)	|	![](sample_output/strain_map/map4.png)