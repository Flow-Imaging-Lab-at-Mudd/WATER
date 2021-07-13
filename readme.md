# Library for manipulating, visualizing, and analyzing synthetic noise of 3D PIV velocity fields

## Data Structure

```matlab
vf = VelocityField(X, U)
```

Is a data structure for a 3D velocity vector field with two essential data/variables: the positions on which velocity measurements were made, `X`, and the corresponding velocity vectors, `U`. Each is a 4D matrix whose first three indices are the grid indices, which specify the position, and in the fourth dimension, the 3-position or 3-velocity vector is stored.

Via various import functions, common formats of data from 3D PIV experiment can be converted to a `VelocityField` object. For example, if the positions and velocities are recorded in components in separate arrays organized on the 3D grid representing the interrogation volume

```matlab
vf = VelocityField.import_grid_separate(xw, yw, zw, uw, vw, zw)
```

This is one of the static methods which adapt data of different original formats into the current object representation.

We may wish to inspect and perform computations on a restricted region of our overall volume, say where the velocity measurements are more reliable. Given a rectangular region specified by the beginning and ending index of each dimension

```matlab
vf.setRange([i_0 i_f; j_0 j_f; k_0 k_f])
```

Such indices may be obtained from a set of mapping functions from position to index such as `vf.getIndex_x()`. Currently, the `vf.getIndex()` methods interpolate to the nearest grid position and wrap around the boundaries, so that even an external position will be assigned a valid index on the grid. To exclude external points, one may use the `vf.inBounds([x y z]')` handle to check for inclusion. Index is the preferred means of access by most methods.

Nonetheless, we may also specify the effective region by position, where the input is a typical range array of dimensions 3 x 2, speicifying the positions of the end points in each dimension. It is not required that the position is given in ascending order, for some fields has descending positions with increasing indices. The indexing is handled properly as to proceed in the same direction of position as in the original field. **It is required however that the subsetted region, and the velocity field in general, be a 3D region--with at least two distinct positions per dimension.**

```matlab
vf.setRangePosition([-5 5; 5 -5; 0 3])
```

## Graphing

We may plot a vector field over our earlier specified range, here the velocity

```matlab
plt = vf.plotVector(vf.U, noise = 0, title_str = 'Velocity $\vec{u}$')
```

Which produces
![global velocity](https://github.com/epicderek/flow/blob/master/illu/3dv.jpg)

The vector field plotted here, `vf.U`, or involved in computation in other methods, if given in the global range, not restricted to the region of interest specified, is automatically subsetted. 

To display the velocity only upon an arbitrary plane perpendicular to the x, y, or z unit vectors, we specify the plane with a normal vector and a base position. Suppose we'd like a plane orthogonal to the x-axis, and this plane is the i<sup>th</sup> such orthogonal plane in on the grid,

```matlab
vf.plotPlaneVector(vf.U, eq = vf.getRegPlaneEq([i 0 0]), noise = 0, title_str = "Velocity $\vec{u}$")
```

Where `vf.getRegPlaneEq([i 0 0])` calculates the required format of the plane as two vectors.

![plane velocity](https://github.com/epicderek/flow/blob/master/illu/plane.jpg)

Now, we introduce noise to our system. This noise will be added to the instance variable `vf.N`, not directly added to `vf.U`, though the user may perform such an addition. Adding Gaussian white noise,

```matlab
vf.noise_wgn(sd = 0, snr = 10)
```

We can visualize the magnitude of this noise along a regular plane, a plane orthogonal to one of the basis vectors. Here we use an additional parameter to specify the plane, considering the planar region probably differs from the 3D region of interest specified earlier. Picking the i<sup>th</sup> orthogonal plane to the x axis, we use `range = [i i; j_0 j_f; k_0 k_f]`, where the dimension with identical beginning and ending indices indicates the direction of normality as well as the index of the plane, and the other two dimensions specify the range of the plane on which noise is to be plotted.

```matlab
vf.plotPlaneScalar(sqrt(sum(vf.N.^2, 4)), range, noise = 0, title_str = 'noise $\delta u$')
```
![plane velocity](https://github.com/epicderek/flow/blob/master/illu/noise_plane.jpg)

To show noise on multiple planes, we generate a Matlab slice plot. Obtaining the equations for the planes in the standard format, stored as a n x 3 x 2 matrix `eqs`, we call

```matlab
vf.slicePlanes(sqrt(sum(vf.N.^2,4)), eqs, noise = 0, 'noise $\Delta u$');
```

![slices](https://github.com/epicderek/flow/blob/master/illu/noise_slice.jpg)

The blank grids are where the planar positions are not proximate enough to the positions on the grid.

Finally, though not applicable to large regions and not especially insightful, we can make a scatter plot of the noise.

```matlab
 vf.plotScalar(sqrt(sum(vf.N.^2, 4)), noise = '$\Delta u$')
```
![scatter3](https://github.com/epicderek/flow/blob/master/illu/scalar-scatter.jpg)

For a continuous scalar field, in addition to `plotPlaneScalar()` and `slicePlanes()`, a standard isosurface plot can also be generated. The required proximity of the actual values in the field to the values specified is as ordained by Matlab's `isosurface()`, which may be examined further. Here we plot two isosurfaces of speed.

```matlab
vf.isosurfaces(vf.data.speed, [250, 200], 0, '$u$')
```

![isosurface](https://github.com/epicderek/flow/blob/master/illu/isosurface.jpg)

## Integration on a Cubic Surface

Surface integration is supported on a rectangular surface, specified by setting the effective region. Using as an example a synthetic Hill's vortex, a spherical structure, we compute the mass flux on the cubic surface `[-1 1; -1 1; -1 1]` which envelopes it.



