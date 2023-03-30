#ifndef __THINLENS__
#define __THINLENS__


#include <vector>
#include <iostream>
#include <fstream>
#include "Camera.h"
#include "GeometricObject.h"
#include "Points.h"
#include "Vectors.h"
#include "Light.h"
#include "Ambient.h"
#include "World.h"


class ThinLensCamera: public Camera {
    public:
        ThinLensCamera(int _hr, int _vr, double d, const Vec3D& _up, const Point3D& pos, const Point3D &_lookAt, float p, int s, double f) 
        {
            hr = _hr;
            vr = _vr;
            distance = d;
            up = _up;
            cameraPos = pos;
            lookAt = _lookAt;
            pixelSize = p;
            samples = s;
            focalPlaneDistance = f;
        }
        ~ThinLensCamera() {}

    private:
        void set_sampler();
        void render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient);

};


void ThinLensCamera::render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient)
{
    Vec3D toPixel = w*distance + right*(-pixelQtnH/2.0) + iup*(pixelQtnV/2.0);/* - (camera.iup/2.0) + (camera.right/2.0)*/ //while using anti-aliasing there is no need to be in the center of the pixel
    Point3D screenP = cameraPos + toPixel;
    Vec3D down;
    std::vector<Vec3D> pixels;
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++)
    {
        if ((i) % (int)pixelQtnH == 0)
        {
            down = down - iup;
            screenP = cameraPos + toPixel;
            screenP = screenP + down;
        } else {
            screenP = screenP + right;
        }
        //anti-aliasing
        int samplesByRow = sqrt(sampler_ptr->get_num_samples());
        Vec3D sum;
        for(int iSamples = 0; iSamples < sampler_ptr->get_num_samples(); iSamples++)
        {  
            Point2D aliasUnit = sampler_ptr->sample_unit_square();
            Vec3D sampleX = right*aliasUnit.x;
            Vec3D sampleY = iup*(-1)*(aliasUnit.y);
            Point3D aliasP = screenP + sampleX + sampleY;

            Point2D diskUnit = sampler_ptr->sample_unit_disk();
            Point3D lensSample = cameraPos + u*((diskUnit.x*2.0) - 1) + v*((diskUnit.y*2.0) - 1);

            double px = aliasP.x * (focalPlaneDistance/distance);
            double py = aliasP.y * (focalPlaneDistance/distance);
            Vec3D focalDir = (w*focalPlaneDistance) + (v*py) + (u*px);
            Point3D onFocalPlane = cameraPos + focalDir;
            sum = sum + trace(lensSample, onFocalPlane, objetos, (*this), lights, &ambient, ambient.depth);    
        }
        pixels.push_back(sum/(double)(sampler_ptr->get_num_samples()));
    }
    std::ofstream pixelOutput("./image.ppm", std::ios::out | std::ios::binary);
    pixelOutput << "P6\n" << pixelQtnH << " " << pixelQtnV << "\n255\n";
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++)
    {
        pixelOutput <<(unsigned char)(std::max(double(1), pixels[i].x)) <<
            (unsigned char)(std::max(double(1), pixels[i].y)) <<
            (unsigned char)(std::max(double(1), pixels[i].z));
    }
    pixelOutput.close();
}


void ThinLensCamera::set_sampler()
{
    sampler_ptr = new JitteredSampler(samples);
}



#endif