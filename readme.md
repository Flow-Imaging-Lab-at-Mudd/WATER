Library for manipulating and visualizing 3D PIV velocity fields.

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

Such indices may be obtained from a set of mapping functions from position to index such as `vf.getIndex_x()`. Indices is the preferred means of access by most methods.

We may plot a vector field over our earlier specified range, here the velocity

```matlab
plt = vf.plotVector(vf.U, noise = 0, title_str = 'Velocity $\vec{u}$')
```

Producing
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
vf.plotPlaneScalar(sqrt(sum(vf.N.^2, 4)), noise = 0, range, title_str = 'noise $\delta u$')
```

(Figure goes here)

To show noise on multiple planes, we generate a Matlab slice plot.

