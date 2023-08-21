#ifndef __PINHOLE__
#define __PINHOLE__


#include <vector>
#include "Camera.h"
#include "JitteredSampler.h"
#include "World.h"
#include "Points.h"
#include "Vectors.h"


class PinholeCamera: public Camera {
    public:
        PinholeCamera(int _hr, int _vr, double d, const Vec3D& _up, const Point3D& pos, const Point3D &_lookAt, float p, int s) {
            hr = _hr;
            vr = _vr;
            distance = d;
            up = _up;
            cameraPos = pos;
            lookAt = _lookAt;
            pixelSize = p;
            samples = s;
        }
        ~PinholeCamera() {}
    private:
        void set_sampler();
        void render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient);
};

void PinholeCamera::render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient)
{
    Vec3D toPixel = w*distance + right*(-pixelQtnH/2.0) + iup*(pixelQtnV/2.0);/* - (camera.iup/2.0) + (camera.right/2.0)*/ //while using anti-aliasing there is no need to be in the center of the pixel
    Vec3D down;
    Vec3D dir;
    std::vector<RGBColor> pixels;
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++)
    {
        if ((i) % (int)pixelQtnH == 0)
        {
            down = down - iup;
            dir = toPixel + down;
        } else {
            dir = dir + right;
        }
        if (i == 19900) {
            int pou = 0;
        }
        //anti-aliasing
        int samplesByRow = sqrt(sampler_ptr->get_num_samples());
        RGBColor sum;
        for(int iSamples = 0; iSamples < sampler_ptr->get_num_samples(); iSamples++)
        {  
            Point2D aliasUnit = sampler_ptr->sample_unit_square();
            Vec3D sampleX = right*aliasUnit.x;
            Vec3D sampleY = iup*(-1)*(aliasUnit.y);
            Vec3D aliasDir = dir + sampleX + sampleY;
            sum = sum + trace(Ray(cameraPos, aliasDir), objetos, (*this), lights, &ambient, ambient.depth);    
        }
        pixels.push_back(sum/(double)(sampler_ptr->get_num_samples()));
    }
    std::ofstream pixelOutput("./image.ppm", std::ios::out | std::ios::binary);
    pixelOutput << "P6\n" << pixelQtnH << " " << pixelQtnV << "\n255\n";
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++)
    {
        pixelOutput <<(unsigned char)(std::max(double(1), pixels[i].r)) <<
            (unsigned char)(std::max(double(1), pixels[i].g)) <<
            (unsigned char)(std::max(double(1), pixels[i].b));
    }
    pixelOutput.close();
}

void PinholeCamera::set_sampler()
{
    sampler_ptr = new JitteredSampler(samples);
}

RGBColor PinholeCamera::setBackgroundSmoothness(const Point3D &pixel, Camera *camera) 
{
    Point2D screenCoordinates = camera->worldToScreenCoordinates(pixel, this->cameraPos);
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




#endif