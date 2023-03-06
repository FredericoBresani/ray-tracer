#ifndef __SAMPLER__
#define __SAMPLER__

#include "Points.h"
#include <vector>

class Sampler 
{
    public:
        virtual void generate_samples(void) = 0;
        void setup_shuffled_indices(void);
        int get_num_samples(void);
        void shuffle_samples(void);
        Point2D sample_unit_square(void);

    protected:
        int n_samples = 1;
        int n_sets = 1;
        std::vector<Point2D> samples;
        std::vector<int> shuffled_indices;
        unsigned long count = 0;
        int jump;
};

Point2D Sampler::sample_unit_square(void)
{
    if (count % n_samples == 0)
    {
        jump = ((int)rand() % n_sets) * n_samples;
    }
    return (samples[(jump + count++) % (n_samples * n_sets)]);
}

int Sampler::get_num_samples()
{
    return n_samples;
}

#endif