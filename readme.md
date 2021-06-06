Library for manipulating and visualizing 3D PIV velocity fields.

```matlab
vf = VelocityField(X, U)
```

Is a data structure for a 3D velocity vector field with two essential data/variables: the positions on which velocity measurements were made, X, and the corresponding velocity vectors, U. Each is a 4D matrix whose first three indices are the grid indices, which specify the position, and in the fourth dimension, the 3-position or 3-velocity vector is stored.

Via various import functions, common formats of data from 3D PIV experiment can be converted to a `VelocityField` object. For example,

```matlab
vf = import_grid_separate(xw, yw, zw, uw, vw, zw)
```

The velocity field can be visualized over the entire grid, adding no noise, as

```matlab
plt = vf.plotVelocity(with_noise = false)
```

Which produces
![global velocity](https://github.com/epicderek/flow/blob/master/illu/3dv.jpg)

To display the velocity only along a plane perpendicular to the x, y, or z unit vectors, we specify the plane with a 3-vector. Suppose we'd like a plane orthogonal to the x-axis, and this plane is the i<sup>th</sup> such orthogonal plane in on the grid,

```matlab
vf.plotPlane(vf.U, index = [i 0 0], noise = 0, title_str = "Velocity $\vec{u}$")
```

Which subsets the velocity field and shows only the velocities on the plane.
![plane velocity](https://github.com/epicderek/flow/blob/master/illu/plane.jpg)

We can also plot the velocity vectors along an arbitrary plane given a unit normal vector and a base position. For a normal vector `orth` and a position `base` stored as column vectors in a matrix `eq`, we use

```matlab
vf.plotPlaneSkewed(vf.U, x = [], eq, noise = 0, title_str = "Velocity $\vec{u}$")
```
(the `x` parameter also allows for specifying three non-colinear points on the plane as columns in a matrix)

Which shows 
![skewed plane velocity](https://github.com/epicderek/flow/blob/master/illu/plane_skew.jpg)


