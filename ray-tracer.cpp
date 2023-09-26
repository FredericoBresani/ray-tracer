#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <stdlib.h>

#include "./header-files/Definitions.h"
#include "./header-files/Matrices.h"
#include "./header-files/Vectors.h"
#include "./header-files/Points.h"
#include "./header-files/RGBColor.h"
#include "./header-files/Ray.h"
#include "./header-files/GeometricObject.h"
#include "./header-files/Sphere.h"
#include "./header-files/Triangle.h"
#include "./header-files/TriangleMesh.h"
#include "./header-files/Plane.h"
#include "./header-files/Line.h"
#include "./header-files/Ambient.h"
#include "./header-files/Light.h"
#include "./header-files/Camera.h"
#include "./header-files/PinholeCamera.h"
#include "./header-files/ThinLensCamera.h"
#include "./header-files/FishEyeCamera.h"
#include "./header-files/PointLight.h"
#include "./header-files/HitInfo.h"
#include "./header-files/Material.h"


// lets use doubles for object-ray intersection and floats for shading calculations

RGBColor setPixelColorNormal(const Vec3D &normal)
{
    return RGBColor(std::min(double(1), normal.x), std::min(double(1), normal.y), std::min(double(1), -normal.z))*210.0;
}

RGBColor setPixelColorCoordinates(const Point3D &location)
{   
    //x = red, y = green, z = blue
    double red = 0.0, green = 0.0, blue = 0.0;
    RGBColor aux = RGBColor(location.x, location.y, location.z)/255.0;
    red += aux.r < 0.0 ? 0.0 : aux.r;
    green += aux.g < 0.0 ? 0.0 : aux.g;
    blue += aux.b < 0.0 ? 0.0 : aux.b;
    aux = RGBColor(std::min(double(1), red*40.0), std::min(double(1), green*40.0), std::min(double(1), blue*40.0))*255.0;
    return aux;
}

/*

RGBColor setBackgroundRGBCoordinates(const Point3D &pixel, Camera *camera) 
{
    Point2D screenCoordinates = camera->worldToScreenCoordinates(pixel, camera->cameraPos);
    double maxScreen = 0, x = 0, y = 0;
    if (camera->hr >= camera->vr) {
        maxScreen = double(camera->hr)/double(camera->vr);
        x = (maxScreen + screenCoordinates.x)/(2.0*maxScreen);
        y = (1.0 + screenCoordinates.y)/2.0;
    } else {
        maxScreen = double(camera->vr)/double(camera->hr);
        x = (1.0 + screenCoordinates.x)/2.0;
        y = (maxScreen + screenCoordinates.y)/(2.0*maxScreen);
    }
    
    return RGBColor(255.0*x, 255.0*y, 60.0);
}*/

void render(std::vector<Object*>& objetos, std::vector<Light*>& lights, Camera& camera, Ambient& ambient)
{
    camera.render(objetos, lights, ambient);
}

int main()
{
    Matrix4D matrix = Matrix4D();
    std::vector<Object*> objetos;
    std::vector<Light*> lights;
    Camera *camera;
    Ambient *ambient;
    char objectType;
    bool read = true;
    float _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18;
    while (scanf("%c %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &objectType, &_1, &_2, &_3, &_4, &_5, &_6, &_7, &_8, &_9, &_10, &_11, &_12, &_13, &_14, &_15, &_16, &_17, &_18) != EOF && read) 
    {
        switch(objectType)
        {
            case 's': 
            {
                Material *mater = new Material{
                    RGBColor(_5, _6, _7), _8, _9, _10, _11, _12, _13, (bool)_15, _16
                };
                Sphere *e = new Sphere(Point3D(_1, _2, _3), _4, mater, (bool)_14);
                if (objetos.size() == 0) {
                    matrix.matrix[0] = {1.0, 0.0, 0.0, 1.0};
                    matrix.matrix[1] = {0.0, 1.0, 0.0, 1.0};
                    matrix.matrix[2] = {0.0, 0.0, 1.0, 0.0};
                    matrix.matrix[3] = {0.0, 0.0, 0.0, 1.0};
                    //e->center = e->center*matrix;
                }
                objetos.push_back(e);
                break;    
            }
            case 'p':
            {
                Material *mater = new Material{
                    RGBColor(_7, _8, _9), _10, _11, _12, _13, _14, _15, (bool)_17, _18
                };
                Plane *p = new Plane(Vec3D(_4, _5, _6), Point3D(_1, _2, _3), mater, (bool)_16);
                objetos.push_back(p);
                break;
            }
            case 't':
            {
                Material *mater = new Material{
                    RGBColor(_3, _4, _5), _6, _7, _8, _9, _10, _11, (bool)_13, _14
                };
                TriangleMesh *mesh = new TriangleMesh((int)_1, (int)_2, mater, (bool)_12);

                float v1, v2, v3;
                int i1, i2, i3;
                while(_2 > 0)
                {
                    scanf("%f %f %f\n", &v1, &v2, &v3);
                    mesh->vertices.push_back(Point3D(v1, v2, v3));
                    _2--;
                }

                while(_1 > 0)
                {
                    scanf("%i %i %i\n", &i1, &i2, &i3);
                    mesh->triangles.push_back(Point3I(i1, i2, i3));
                    _1--;
                }
                
                objetos.push_back(mesh);
                break;
            }
            case 'l':
            {
                PointLight *l = new PointLight(Point3D(_1, _2, _3), RGBColor(_4, _5, _6), (bool)_7);
                lights.push_back(l);
                break;
            }
            case 'c':
            {
                camera = new PinholeCamera(_1, _2, _3, Vec3D(_4, _5, _6), Point3D(_7, _8, _9), Point3D(_10, _11, _12), _13, _14);
                // camera = new ThinLensCamera(_1, _2, _3, Vec3D(_4, _5, _6), Point3D(_7, _8, _9), Point3D(_10, _11, _12), _13, _14, _15);
                // camera = new FishEyeCamera(_1, _2, _3, Vec3D(_4, _5, _6), Point3D(_7, _8, _9), Point3D(_10, _11, _12), _13, _14, _15);


                //1rad = 180/pi graus
                double cos = std::cos(M_PI/6.0); 
                double sen = std::sin(M_PI/6.0);
                matrix.matrix[0] = {cos, 0.0, -sen, 1.0};
                matrix.matrix[1] = {0.0, 1.0, 0.0, 2.0};
                matrix.matrix[2] = {sen, 0.0, cos, -1.0};
                matrix.matrix[3] = {0.0, 0.0, 0.0, 1.0};
                //camera->transformCamera(matrix);
                camera->makeCamera();
                break;
            }
            case 'r':
            {
                Line *line = new Line(Point3D(_1, _2, _3), Vec3D(_4, _5, _6), RGBColor(_7, _8, _9));
                // line->setNormal(camera->cameraPos);
                objetos.push_back(line);
                break;
            }  
            case 'a':
            {
                ambient = new Ambient(RGBColor(_1, _2, _3), _4, _5);
                break;
            }
            default:
            {
                read = false;
                break;
            }
        }
    }
    render(objetos, lights, *camera, *ambient);
    return 0;
}