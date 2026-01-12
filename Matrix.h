#pragma once
#include <cuda_runtime.h>

template <class T>
class Matrix {
private:
    int rows;
    int cols;
    T* data;
public:
    Matrix(int m_rows, int m_cols) : rows(m_rows), cols(m_cols) { cudaMallocManaged(&data, m_rows * m_cols * sizeof(T)); }
    ~Matrix() { cudaFree(data); }
    T& operator()(int i) { return data[i]; }
    T& operator()(int i, int j) { return data[i * cols + j]; }
    T* data_ptr() { return data; }
};