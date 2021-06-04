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
plt = vf.plotVelocity(noise = 0)
```

Which produces
![Image of Yaktocat](https://octodex.github.com/images/yaktocat.png)

