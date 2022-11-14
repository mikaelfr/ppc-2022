#include <cstdlib>
#include <cmath>
#include <iostream>

// cores: 10 and 20 for remote, 4 and 8 for local
const int NUM_PARALLEL = 8;
const int NUM_INST = 4;

#define SET_INST(name) double name ## _ [NUM_INST] = {};
#define RUN_INST(name, x) for (int a = 0; a < NUM_INST; a++) { x } 
#define AGG_INST(name) for (int a = 0; a < NUM_INST; a++) { name += name ## _ [a]; }

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
    int na = nx / NUM_INST;
    int nab = na * NUM_INST;

    double* d = (double*)std::malloc(nx * ny * sizeof(double));

    for (int y = 0; y < ny; y++) {
        double sum = 0.0;
        SET_INST(sum);
        for (int x = 0; x < nab; x += NUM_INST) {
            RUN_INST(sum,
                sum_[a] += data[x+a + y*nx];
            )
        }
        for (int x = nab; x < nx; x++){
            sum += data[x + y*nx];
        }
        AGG_INST(sum);
        sum /= nx;

        double sqsum = 0.0;
        SET_INST(sqsum);
        for (int x = 0; x < nab; x += NUM_INST) {
            RUN_INST(sqsum,
                double val = (double)data[x+a + y*nx] - sum;
                d[x+a + y*nx] = val;
                sqsum_[a] += val * val;
            )
        }
        for (int x = nab; x < nx; x++){
            double val = (double)data[x + y*nx] - sum;
            d[x + y*nx] = val;
            sqsum += val * val;
        }
        AGG_INST(sqsum);
        sqsum = std::sqrt(sqsum);

        for (int x = 0; x < nx; x++) {
            d[x + y*nx] /= sqsum;
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double cor = 0.0;
            SET_INST(cor);
            for (int x = 0; x < nab; x += NUM_INST) {
                RUN_INST(cor,
                    cor_[a] += d[x+a + i*nx] * d[x+a + j*nx];
                )
            }
            for (int x = nab; x < nx; x++){
                cor += d[x + i*nx] * d[x + j*nx];
            }
            AGG_INST(cor);
            result[i + j*ny] = cor;
        }
    }

    std::free(d);
}
