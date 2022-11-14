#include <cstdlib>
#include <cmath>
#include <iostream>

const int DVSIZE = 8;

// The grading machine has 512bit wide vectors
typedef double double8_t __attribute__ ((vector_size (8 * sizeof(double))));

constexpr double8_t d8zero {
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
};

#define SET_VEC(name) double8_t name ## 8 = d8zero;
#define RUN_VEC(name, x) for (int v = 0; v < DVSIZE; v++) { x } 
#define AGG_VEC(name) for (int v = 0; v < DVSIZE; v++) { name += name ## 8 [v]; }

static double8_t* double8_alloc(std::size_t n) {
    void* tmp = 0;
    if (posix_memalign(&tmp, sizeof(double8_t), sizeof(double8_t) * n)) {
        throw std::bad_alloc();
    }
    return (double8_t*)tmp;
}

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
    int na = (nx + DVSIZE - 1) / DVSIZE;
    int nab = na * DVSIZE;
    int remi = DVSIZE - (nab - nx);

    double8_t* vd = double8_alloc(ny*na);

    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < na; x++) {
            for (int i = 0; i < DVSIZE; i++) {
                const int ox = x * DVSIZE + i;
                vd[x + y*na][i] = (ox < nx) ? data[ox + y*nx] : 0.0;
            }
        }

        double sum = 0.0;
        double8_t sum8 = d8zero;
        for (int x = 0; x < na; x++) {
            sum8 += vd[x + y*na];
        }
        AGG_VEC(sum);
        sum /= nx;

        // filling in the couple last elements of the last vector with the sum
        // so when we normalize, these go properly to 0
        for (int vx = remi; vx < DVSIZE; vx++) {
            vd[na-1 + y*na][vx] = sum;
        }

        double sqsum = 0.0;
        double8_t sqsum8 = d8zero;
        for (int x = 0; x < na; x++) {
            double8_t val = vd[x + y*na] - sum;
            vd[x + y*na] = val;
            sqsum8 += val * val;
            
        }
        AGG_VEC(sqsum);
        sqsum = std::sqrt(sqsum);

        for (int x = 0; x < na; x++) {
            vd[x + y*na] /= sqsum;
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double cor = 0.0;
            double8_t cor8 = d8zero;
            for (int x = 0; x < na; x++) {
                cor8 += vd[x + i*na] * vd[x + j*na];
            }
            AGG_VEC(cor);
            result[i + j*ny] = cor;
        }
    }

    std::free(vd);
}
