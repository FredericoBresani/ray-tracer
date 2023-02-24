# ray-tracer

#### input format
- spheres:
```bash
s c1 c2 c3 r R G B kd ks ka kr kt p
# s: identify a sphere
# (c1, c2, c3): the coordinates of the sphere center
# r: sphere radius
# (R, G, B): the sphere difuse color
```
- planes:
```bash
p x y z v1 v2 v3 R G B kd ks ka kr kt p
# p: identify a plane
# (x, y, z): the coordinates of a point in the plane
# (v1, v2, v3): the coordinates of the plane normal
```

- triangles:
```bash
t a1 a2 a3 b1 b2 b3 c1 c2 c3 R G B kd ks ka kr kt p
# t: identify a triangle
# (a1, a2, a3): the coordinates of the first vertex of the triangle
# (b1, b2, b3): the coordinates of the second vertex of the triangle
# (c1, c2, c3): the coordinates of the third vertex of the triangle
```

- all previous objects:
```bash
# all the objects has the following attributes to define their material
# kd: Difuse coeficient
# ks: Specular coeficient
# ka: Ambient coeficient
# kr: Reflective coeficient
# kt: Transmission coeficient
# p:  Phong exponent
```

- lights:
```bash
l l1 l2 l3 R G B
# l: identify a light
# (l1, l2, l3): the coordinates of the light location
# (R, G, B): the light intensity
```

- camera:
```bash
c h_res v_res d up1 up2 up3 l1 l2 l3 m1 m2 m3
# c: identify a camera
# h_res: horizontal resolution
# v_res: vertical resolution
# d: distance to the screen
# (up1, up2, up3): the coordinates of the up vector
# (l1, l2, l3): the coordinates of the camera location
# (m1, m2, m3): the coordintes of the location the the camera points at
```