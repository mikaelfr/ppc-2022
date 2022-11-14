#include <cstdlib>
#include <cmath>
#include <iostream>

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
void correlate(int ny, int nx, const float *data, float *result) {
    double* d = (double*)std::malloc(ny * nx * sizeof(double));
 
    for (int y = 0; y < ny; y++) {
        double sum = 0.0;
        for (int x = 0; x < nx; x++) {
            sum += data[x + y*nx];
        }
        sum /= nx;

        double sqsum = 0.0;
        for (int x = 0; x < nx; x++) {
            double val = (double)data[x + y*nx] - sum;
            d[x + y*nx] = val;
            sqsum += val * val;
        }
        sqsum = std::sqrt(sqsum);

        for (int x = 0; x < nx; x++) {
            d[x + y*nx] /= sqsum;
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double cor = 0.0;
            for (int x = 0; x < nx; x++) {
                cor += d[x + i*nx] * d[x + j*nx];
            }
            result[i + j*ny] = cor;
        }
    }

    std::free(d);
}
