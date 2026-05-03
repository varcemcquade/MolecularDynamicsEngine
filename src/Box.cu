#include <cuda_runtime.h>

__global__
void ComputeForcesKernel(double* device_coordinates, double* device_forces, const int ATOM_COUNT,
                         const double EPSILON_ANGSTROMS, const double SIGMA,
                         const double BOX_LENGTH, const double CUTOFF_RANGE,
                         const bool isPeriodic) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= ATOM_COUNT) { return; }

    device_forces[i * 3] = 0;
    device_forces[i * 3 + 1] = 0;
    device_forces[i * 3 + 2] = 0;

    double xi = device_coordinates[i * 3];
    double yi = device_coordinates[i * 3 + 1];
    double zi = device_coordinates[i * 3 + 2];

    double SIGMA3 = SIGMA * SIGMA * SIGMA;
    double SIGMA6 = SIGMA3 * SIGMA3;
    double SIGMA12 = SIGMA6 * SIGMA6;

    for (int j = 0; j < ATOM_COUNT; ++j) {
        if (i == j) { continue; }

        double dx = xi - device_coordinates[j * 3];
        double dy = yi - device_coordinates[j * 3 + 1];
        double dz = zi - device_coordinates[j * 3 + 2];

        if (isPeriodic) {
            dx -= BOX_LENGTH * round(dx / BOX_LENGTH);
            dy -= BOX_LENGTH * round(dy / BOX_LENGTH);
            dz -= BOX_LENGTH * round(dz / BOX_LENGTH);
        }

        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);

        if (r > CUTOFF_RANGE) { continue; }

        double r6 = r2 * r2 * r2;
        double r8 = r6 * r2;
        double r14 = r8 * r6;

        double u = ((6 * SIGMA6) / r8) - ((12 * SIGMA12) / r14);
        double scalar_f = -4 * EPSILON_ANGSTROMS * u;

        device_forces[i * 3] += scalar_f * dx;
        device_forces[i * 3 + 1] += scalar_f * dy;
        device_forces[i * 3 + 2] += scalar_f * dz;
    }
}

extern "C" void AllocateGPUBuffers(double** device_coordinates, double** device_forces,
                                    const int ATOM_COUNT) {
    int size_in_bytes = 3 * ATOM_COUNT * sizeof(double);
    cudaMalloc(reinterpret_cast<void**>(device_coordinates), size_in_bytes);
    cudaMalloc(reinterpret_cast<void**>(device_forces),      size_in_bytes);
}

extern "C" void FreeGPUBuffers(double* device_coordinates, double* device_forces) {
    cudaFree(device_coordinates);
    cudaFree(device_forces);
}

extern "C" void UpdateForcesGPU(double* device_coordinates, double* device_forces,
                                  double* host_coordinates,   double* host_forces,
                                  const int ATOM_COUNT,
                                  const double EPSILON_ANGSTROMS, const double SIGMA,
                                  const double BOX_LENGTH, const double CUTOFF_RANGE,
                                  const bool isPeriodic) {
    int threads       = 256;
    int blocks        = (ATOM_COUNT + threads - 1) / threads;
    int size_in_bytes = 3 * ATOM_COUNT * sizeof(double);

    cudaMemcpy(device_coordinates, host_coordinates, size_in_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(device_forces, host_forces, size_in_bytes, cudaMemcpyHostToDevice);

    ComputeForcesKernel<<<blocks, threads>>>(device_coordinates, device_forces, ATOM_COUNT,
                                             EPSILON_ANGSTROMS, SIGMA,
                                             BOX_LENGTH, CUTOFF_RANGE,
                                             isPeriodic);

    cudaDeviceSynchronize();
    cudaMemcpy(host_coordinates, device_coordinates, size_in_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_forces, device_forces, size_in_bytes, cudaMemcpyDeviceToHost);
}