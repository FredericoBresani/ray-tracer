#ifndef __JITTERED__
#define __JITTEREF__

#include "Sampler.h"
#include "Points.h"
#include <stdlib.h>

class JitteredSampler: public Sampler
{
    public:
        JitteredSampler(const int n)
        {
            n_samples = n*n;
            n_sets = 83;
            this->generate_samples();
        }
        ~JitteredSampler();
    private:
        void generate_samples();
};

void JitteredSampler::generate_samples(void)
{
    int s = sqrt(n_samples);
    for (int sets = 0; sets < n_sets; sets++)
        for (int rows = 0; rows < s; rows++)
            for (int columns = 0; columns < s; columns++)
            {
                double randX = ((double)(columns + rand()))/(double)(RAND_MAX + columns);
                double randY = ((double)(rows + rand()))/(double)(RAND_MAX + rows); 
                samples.push_back(Point2D(randX, randY));
            }
}


#endif