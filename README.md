# ray-tracer

### Some rendered examples

#### Colored shadows
![Colored shadows](presets/colored-shadows.png)
- This example has no kinds of reflection;
-----

#### Spheres in limbo
![Spheres in limbo](presets/spheres-in-limbo.png)
- Here I get some diffuse reflection;
-----
#### Transparent, metalic, glass and plastic spheres
![Transparent, metalic, glass and plastic sphere](presets/spheres-transp-metal-glass.png)
- No diffuse reflection;
- Reflections and transparency;
-----
#### Thins lens camera
![Thin lens camera](presets/thin-lens-camera.png)
- With this implementation of the camera I'm able to render imagens choosing focal planes;
- Here I'm using triangle meshes;
-----
#### Metalic spheres and triangles inside walls
![Metalic spheres and triangles inside walls](presets/spheres-triangles-walls.png)
- Here I got a noisy image do to the diffuse reflection;

-----

#### Fish eye camera
![Fish eye camera](presets/fish-eye-example.png)
- A simple implementation of a fish eye camera;
### input format
#### Spheres
```bash
s c1 c2 c3 r R G B kd ks ka kr kt p s ms ior
# s: identify a sphere
# (c1, c2, c3): the coordinates of the sphere center
# r: sphere radius
# (R, G, B): the sphere difuse color
```

-----

#### Planes
```bash
p x y z v1 v2 v3 R G B kd ks ka kr kt p s ms ior
# p: identify a plane
# (x, y, z): the coordinates of a point in the plane
# (v1, v2, v3): the coordinates of the plane normal
```

-----


#### Triangle mesh
```bash
t nt nv R G B kd ks ka kr kt p s ms ior
p1x p1y p1z
p2x p2y p2z
.
.
.
pnvx pnvy pnvz
t1a t1b t1c
t2a t2b t2c
.
.
.
tnta tntb tntc
# t: identify a triangle
# (p1x, p1y, p1z) to (pnvx, pnvy, pnvz): the coordinates of the the vertices
# (t1a, t1b, t1c) to (tnta, tntb, tntc): triples of vertices indices
```

-----


#### All previous objects:
```bash
# all the objects has the following attributes to define their material
# kd: Difuse coeficient
# ks: Specular coeficient
# ka: Ambient coeficient
# kr: Reflective coeficient
# kt: Transmission coeficient
# p:  Phong exponent
# s: boolean to tell if the object could cast shadows
# ms: boolean inside the object material to know if the material can get shadowed
# ior: refraction coeficient
```

-----


#### Lights
```bash
l l1 l2 l3 R G B s
# l: identify a light
# (l1, l2, l3): the coordinates of the light location
# (R, G, B): the light intensity
# s: boolean to tell if the light could produce shadows
```

-----

#### Camera
```bash
c h_res v_res d up1 up2 up3 l1 l2 l3 m1 m2 m3 p s optional
# c: identify a camera
# h_res: horizontal resolution
# v_res: vertical resolution
# d: distance to the screen
# (up1, up2, up3): the coordinates of the up vector
# (l1, l2, l3): the coordinates of the camera location
# (m1, m2, m3): the coordintes of the location the the camera points at
# p: pixel size
# s: samples
# optional: it could be the focal plane distance or the fish eye camera max angle (from 0 to 2PI rad)
```

-----

#### Ambient:
```bash
a R G B ir depth
# (R, G, B): The ambient light color
# ir: The ambient reflectiveness coeficient
# depth: This is the depth the ray-tracer will use
```