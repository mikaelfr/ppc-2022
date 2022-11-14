#include <cstdlib>
#include <cmath>
#include <iostream>

const int DVSIZE = 8;
// two 512bit vector ops per core
const int NUM_INST = 4;

// The grading machine has 512bit wide vectors
typedef double double8_t __attribute__ ((vector_size (8 * sizeof(double))));

constexpr double8_t d8zero {
    0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0
};

#define PRINT_VAL(name) std::cout << #name << ": " << name << std::endl;
#define PRINT_ARRAY() for (int x = 0; x < nabi; x++) { for (int i = 0; i < DVSIZE; i++) { std::cout << vd[x + y*nabi][i] << " "; }}; std::cout << std::endl;

#define SET_VEC(name) double8_t name ## 8 = d8zero;
#define RUN_VEC(name, x) for (int v = 0; v < DVSIZE; v++) { x } 
#define AGG_VEC(name) for (int v = 0; v < DVSIZE; v++) { name += name ## 8 [v]; }

#define SET_INST(name) double8_t name ## _ [NUM_INST] = { };
#define RUN_INST(name, x) for (int a = 0; a < NUM_INST; a++) { x } 
#define AGG_INST(name) for (int a = 0; a < NUM_INST; a++) { name += name ## _ [a]; }

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
    const int na = (nx + DVSIZE - 1) / DVSIZE;
    const int nab = na * DVSIZE;
    
    const int nai = (na + NUM_INST - 1) / NUM_INST;
    // number of padded vectors when taking into account 
    // number of vectors (na) has to be divisible by 2
    // vectors are length DVSIZE
    const int nabi = nai * NUM_INST;

    // all remaining values
    const int rem = nabi*DVSIZE - nx;
    // remaining number of vectors (partly filled or empty)
    const int remv = (rem + DVSIZE) / DVSIZE;
    // starting index for the partly filled vector
    const int remvi = nabi - remv;
    // starting index in the partly filled vector
    const int remi = DVSIZE - (nab - nx);

    double8_t* vd = double8_alloc(ny*nabi);

    #pragma omp parallel for schedule(static, 1)
    for (int y = 0; y < ny; y++) {
        for (int x = 0; x < nabi; x++) {
            for (int i = 0; i < DVSIZE; i++) {
                const int ox = x * DVSIZE + i;
                vd[x + y*nabi][i] = (ox < nx) ? data[ox + y*nx] : 0.0;
            }
        }

        double sum = 0.0;
        double8_t sum8 = d8zero;
        SET_INST(sum8);
        for (int x = 0; x < nabi; x += NUM_INST) {
            RUN_INST(sum,
                sum8_[a] += vd[x+a + y*nabi];
            )
        }
        AGG_INST(sum8);
        AGG_VEC(sum);
        sum /= nx;

        // filling in the couple last elements of the last partly filled vector
        // with the sum so when we normalize, these go properly to 0
        for (int vx = remi; vx < DVSIZE; vx++) {
            vd[remvi + y*nabi][vx] = sum;
        }

        // filling in the completely empty vectors with sum
        for (int x = remvi + 1; x < nabi; x++) {
            for (int i = 0; i < DVSIZE; i++) {
                vd[x + y*nabi][i] = sum;
            }
        }

        double sqsum = 0.0;
        double8_t sqsum8 = d8zero;
        SET_INST(sqsum8);
        for (int x = 0; x < nabi; x += NUM_INST) {
            RUN_INST(sqsum8,
                double8_t val = vd[x+a + y*nabi] - sum;
                vd[x+a + y*nabi] = val;
                sqsum8_[a] += val * val;
            )
        }
        AGG_INST(sqsum8);
        AGG_VEC(sqsum);
        sqsum = std::sqrt(sqsum);

        for (int x = 0; x < nabi; x++) {
            vd[x + y*nabi] /= sqsum;
        }
    }

    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double cor = 0.0;
            double8_t cor8 = d8zero;
            SET_INST(cor8);
            for (int x = 0; x < nabi; x += NUM_INST) {
                RUN_INST(cor8,
                    cor8_[a] += vd[x+a + i*nabi] * vd[x+a + j*nabi];
                )
            }
            AGG_INST(cor8);
            AGG_VEC(cor);
            result[i + j*ny] = cor;
        }
    }

    std::free(vd);
}
