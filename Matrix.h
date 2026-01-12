#pragma once
#include <vector>

using namespace std;

template <class T>
class Matrix {
private:
	int rows;
	int cols;
	vector<T> data;
public:
	Matrix();
	Matrix(int m_rows, int m_cols);
	~Matrix();

	T& operator()(int i);
	T& operator()(int i, int j);
};


template <class T>
Matrix<T>::Matrix() : rows(0), cols(0) {};

template <class T>
Matrix<T>::Matrix(int m_rows, int m_cols) : rows(m_rows), cols(m_cols), data(m_rows* m_cols) {};

template <class T>
Matrix<T>::~Matrix() {};

template <class T>
T& Matrix<T>::operator()(int i) { return data[i]; }

template <class T>
T& Matrix<T>::operator()(int i, int j) { return data[i * cols + j]; }