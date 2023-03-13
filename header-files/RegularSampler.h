#ifndef __REGULARSAMPLER__
#define __REGULARSAMPLER__

#include "Sampler.h"
#include "Vectors.h"
#include "Points.h"

class RegularSampler: public Sampler 
{
    public:
        RegularSampler(const int n)
        {
            n_samples = n*n;
            n_sets = 1;
            this->generate_samples();
            this->setup_shuffled_indices();
        }
        ~RegularSampler() {}
        private:
            void generate_samples(void);

};

void RegularSampler::generate_samples(void)
{
    int s = sqrt(n_samples);
    Vec2D right = Vec2D(1.0, 0.0)/(float)s;
    Vec2D down = Vec2D(0.0, 1.0)/(float)s; // the higher the iteration in trace(), the lower down (positive here) gets.
    Point2D firstSample = Point2D() + right/2.0 + down/2.0;
    for(int iSamples = 0; iSamples < s; iSamples++)
    {
        for(int jSamples = 0; jSamples < s; jSamples++)
        {
            samples.push_back(firstSample + right*jSamples + down*iSamples);
        }
    }
}


#endif