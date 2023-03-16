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
#include "./header-files/HitInfo.h"
#include "./header-files/GeometricObject.h"
#include "./header-files/Sphere.h"
#include "./header-files/Triangle.h"
#include "./header-files/Plane.h"
#include "./header-files/Line.h"
#include "./header-files/Ambient.h"
#include "./header-files/Light.h"
#include "./header-files/Camera.h"
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
    RGBColor aux = Vec3D(location.x, location.y, location.z)/255.0;
    red += aux.x < 0.0 ? 0.0 : aux.x;
    green += aux.y < 0.0 ? 0.0 : aux.y;
    blue += aux.z < 0.0 ? 0.0 : aux.z;
    aux = RGBColor(std::min(double(1), red*40.0), std::min(double(1), green*40.0), std::min(double(1), blue*40.0))*255.0;
    return aux;
}

RGBColor setBackgroundSmoothness(const Point3D &pixel, Camera *camera) 
{
    Point2D screenCoordinates = camera->worldToScreenCoordinates(pixel, camera->cameraPos);
    double maxScreen = 0, x = 0, y = 0;
    if (camera->hr >= camera->vr) {
        maxScreen = double(camera->hr)/double(camera->vr);
        y = (1.0 + std::abs(screenCoordinates.y))/2.0;
    } else {
        maxScreen = double(camera->vr)/double(camera->hr);
        y = (maxScreen + std::abs(screenCoordinates.y))/(2.0*maxScreen);
    }
    return RGBColor(135.0, 206.0, 235.0)*y;
}

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
}

RGBColor trace(const Point3D& origin, const Point3D& pixel, std::vector<Object*>& objetos, Camera camera, std::vector<Light*> lights, Ambient *ambient, int depth)
{
    double t = infinity;
    double tmin = infinity;
    double kd, ks, ka, kr, kt, phongExp; 
    HitInfo *hInfo = new HitInfo();
    Ray *ray = new Ray(origin, pixel - origin);
    RGBColor color;
    if (depth == 0) {
        return color;
    }
    for (int i = 0; i < objetos.size(); i++)
    {
        if (objetos[i]->rayObjectIntersect(*ray, &t, *hInfo))
        {
            if (t < tmin)
            {
                tmin = t;
                hInfo->hit_location = ray->origin + ray->direction*tmin;
                hInfo->normal = objetos[i]->getNormal(hInfo->hit_location, *ray);
                color = objetos[i]->getColor();
                kd = objetos[i]->getKd();
                ks = objetos[i]->getKs();
                ka = objetos[i]->getKa();
                kr = objetos[i]->getKr();
                kt = objetos[i]->getKt();
                phongExp = objetos[i]->getPhongExp();
                // color = setPixelColorNormal(hit.normal);
                // color = setPixelColorCoordinates(hit.hit_location);
            }
        }
    }
    
    if (tmin == infinity)
    {
        // return setBackgroundSmoothness(pixel, &camera);
        // return setBackgroundRGBCoordinates(pixel, &camera);
        // return RGBColor(0.0, 0.0, 0.0);
        return RGBColor(255.0, 255.0, 255.0);
        // return RGBColor(190.0, 230.0, 255.0);
    } else {
        double difuseIndice = 0, rMax = 0, gMax = 0, bMax = 0, reflectiveness = 0;
        Vec3D resultingColor, mixedColor;
        RGBColor colorFilter = ambient->color*ka;
        hInfo->toCamera = camera.cameraPos - hInfo->hit_location;
        hInfo->toCamera = hInfo->toCamera.normalize(hInfo->toCamera);
        hInfo->viewerReflex = ((hInfo->normal*2)*(hInfo->normal*hInfo->toCamera)) - hInfo->toCamera;
        hInfo->viewerReflex = hInfo->viewerReflex.normalize(hInfo->viewerReflex);
        Point3D hitPoint = hInfo->hit_location + hInfo->normal*0.001;

        Vec3D auxVec = hInfo->viewerReflex^hInfo->normal;

        Vec3D auxReflex = auxVec^hInfo->viewerReflex;
        Point3D reflexPoint = hitPoint + (hInfo->viewerReflex*10.0) + (auxVec*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + ((auxVec*(-1.0))*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + (auxReflex*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + ((auxReflex*(-1.0))*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr));


        Vec3D auxNormal = auxVec^hInfo->normal;
        Point3D difusePoint = hitPoint + hInfo->normal + (auxVec*((double)std::rand()/(double)RAND_MAX)*kd) + ((auxVec*(-1.0))*((double)std::rand()/(double)RAND_MAX)*kd) + (auxNormal*((double)std::rand()/(double)RAND_MAX)*kd) + ((auxNormal*(-1.0))*((double)std::rand()/(double)RAND_MAX)*kd);
        
        reflectiveness = 1.0 + ambient->ir;
        for (int l = 0; l < lights.size(); l++) {
            hInfo->toLight = lights[l]->lightPos - hInfo->hit_location;
            hInfo->toLight = hInfo->toLight.normalize(hInfo->toLight);
            hInfo->reflection = ((hInfo->normal*2)*(hInfo->normal*hInfo->toLight)) - hInfo->toLight;
            hInfo->reflection = hInfo->reflection.normalize(hInfo->reflection);
            mixedColor = Vec3D(lights[l]->lightColor.x*(color.x*reflectiveness), lights[l]->lightColor.y*(color.y*reflectiveness), lights[l]->lightColor.z*(color.z*reflectiveness))/255.0;
            resultingColor = resultingColor + mixedColor*kd*std::max(hInfo->normal*hInfo->toLight, 0.0);
        }
         
        color = RGBColor(colorFilter.x*resultingColor.x, colorFilter.y*resultingColor.y, colorFilter.z*resultingColor.z)/255.0;
        color = RGBColor(std::min(color.x , 255.0), std::min(color.y, 255.0), std::min(color.z, 255.0));
        if (kr > 0) {
            color = color + trace(hitPoint, reflexPoint, objetos, camera, lights, ambient, depth - 1);
            color = color/2.0;
        } else {
            color = color + trace(hitPoint, difusePoint, objetos, camera, lights, ambient, depth - 1);
            color = color/2.0;
        }
    }
    return color;
}

void render(std::vector<Object*>& objetos, std::vector<Light*>& lights, Camera& camera, Ambient& ambient)
{
    Vec3D toPixel = camera.w*camera.distance + camera.right*(-camera.pixelQtnH/2.0) + camera.iup*(camera.pixelQtnV/2.0);/* - (camera.iup/2.0) + (camera.right/2.0)*/ //while using anti-aliasing there is no need to be in the center of the pixel
    Point3D screenP = camera.cameraPos + toPixel;
    Vec3D down;
    std::vector<Vec3D> pixels;
    for (int i = 0; i < camera.pixelQtnH*camera.pixelQtnV; i++)
    {
        if ((i) % (int)camera.pixelQtnH == 0)
        {
            down = down - camera.iup;
            screenP = camera.cameraPos + toPixel;
            screenP = screenP + down;
        } else {
            screenP = screenP + camera.right;
        }

        //anti-aliasing
        int samplesByRow = sqrt(camera.sampler_ptr->get_num_samples());
        Vec3D sum;
        for(int iSamples = 0; iSamples < camera.sampler_ptr->get_num_samples(); iSamples++)
        {  
            Point2D aliasUnit = camera.sampler_ptr->sample_unit_square();
            Vec3D sampleX = camera.right*aliasUnit.x;
            Vec3D sampleY = camera.iup*(-1)*(aliasUnit.y);
            sum = sum + trace(camera.cameraPos, screenP + sampleX + sampleY, objetos, camera, lights, &ambient, ambient.depth);    
        }

        pixels.push_back(sum/(double)(camera.sampler_ptr->get_num_samples()));
    }
    std::ofstream pixelOutput("./image.ppm", std::ios::out | std::ios::binary);
    pixelOutput << "P6\n" << camera.pixelQtnH << " " << camera.pixelQtnV << "\n255\n";
    for (int i = 0; i < camera.pixelQtnH*camera.pixelQtnV; i++)
    {
        pixelOutput <<(unsigned char)(std::max(double(1), pixels[i].x)) <<
            (unsigned char)(std::max(double(1), pixels[i].y)) <<
            (unsigned char)(std::max(double(1), pixels[i].z));
    }
    pixelOutput.close();   
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
                    RGBColor(_5, _6, _7), _8, _9, _10, _11, _12, _13
                };
                Sphere *e = new Sphere(Point3D(_1, _2, _3), _4, mater);
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
                    RGBColor(_7, _8, _9), _10, _11, _12, _13, _14, _15
                };
                Plane *p = new Plane(Vec3D(_4, _5, _6), Point3D(_1, _2, _3), mater);
                objetos.push_back(p);
                break;
            }
            case 't':
            {
                Material *mater = new Material{
                    RGBColor(_10, _11, _12), _13, _14, _15, _16, _17, _18
                };
                Triangle *t = new Triangle(Point3D(_1, _2, _3), Point3D(_4, _5, _6), Point3D(_7, _8, _9), mater);
                objetos.push_back(t);
                break;
            }
            case 'l':
            {
                Light *l = new Light(Point3D(_1, _2, _3), RGBColor(_4, _5, _6));
                lights.push_back(l);
                break;
            }
            case 'c':
            {
                camera = new Camera(_1, _2, _3, Vec3D(_4, _5, _6), Point3D(_7, _8, _9), Point3D(_10, _11, _12), _13, _14);
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
                line->setNormal(camera->cameraPos);
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