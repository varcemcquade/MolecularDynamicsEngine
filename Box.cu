#ifdef __INTELLISENSE__
#define __global__
#define __device__
#define __host__

struct dim3 {
public:
    int x;
    int y;
    int z;
};

dim3 blockIdx;
dim3 blockDim;
dim3 threadIdx;

#endif

#include <cuda_runtime.h>
#include <cmath>

__global__
void ComputeForcesKernel(double* coords, double* forces, int n, double epsilon, double sigma) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) { return; }

    forces[i * 3] = 0;
    forces[i * 3 + 1] = 0;
    forces[i * 3 + 2] = 0;

    double xi = coords[i * 3];
    double yi = coords[i * 3 + 1];
    double zi = coords[i * 3 + 2];

    for (int j = 0; j < n; ++j) {
        if (i == j) { continue; }

        double dx = xi - coords[j * 3];
        double dy = yi - coords[j * 3 + 1];
        double dz = zi - coords[j * 3 + 2];
        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);

        double sigma3 = sigma * sigma * sigma;
        double sigma6 = sigma3 * sigma3;
        double sigma12 = sigma6 * sigma6;

        double r6 = r2 * r2 * r2;
        double r8 = r6 * r2;
        double r14 = r8 * r6;


        double u = ((6 * sigma6) / r8) - ((12 * sigma12) / r14); // Angstroms^-2
        double scalar_f = -4 * epsilon * u;

        forces[i * 3] += scalar_f * dx;
        forces[i * 3 + 1] += scalar_f * dy;
        forces[i * 3 + 2] += scalar_f * dz;
    }
}

extern "C" void ComputeForcesGPU(double* coordinates, double* forces, int n, double epsilon, double sigma) {
    int threads = 256;
    int blocks = (n + threads - 1) / threads;
    ComputeForcesKernel <<<blocks, threads>>> (coordinates, forces, n, epsilon, sigma);
    cudaDeviceSynchronize();
}