#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/date_time.hpp>
#include <string>
#include <vector>
#include <sstream>
//cul without k sampling


//using namespace boost::property_tree;
//using namespace boost::gregorian;
//using namespace boost;
using namespace std;
//string stra,strc;
//vector<string> vecStr;
//ptree pt,p1,p2;
//stringstream stream;
//const std::string file_path = "./protozoa_k4_t1_v10.json";

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>
#include <queue>//队列先进先出
#include <stack>//栈处理函数
#include <fstream>//定义文件流
#include <sstream>//串流输入输出操作
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <map>
#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"
#include <locale>

#define BASE 4

typedef float REAL;


//---------------------class matrix
#ifndef SIMPLE_MATRIX_H
#define SIMPLE_MATRIX_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <ctime>
#include <cmath>

#define MAX_LINE_LENGTH 100000
#define EPSILON 0.0000001
namespace smat
{

	/**********************************************
	* Declaration part
	**********************************************/

	template<class T>
	class Matrix//定义矩阵类
	{
	public:
		Matrix(int rows, int columns); // initialization without assigning values初始化矩阵不赋值
		Matrix(int rows, int columns, T value); // initialization with all same values初始化矩阵值都为T
		Matrix(int rows, int columns, std::string type); // special matrix such as I特殊矩阵I
		Matrix(const char * filename); // load matrix from txt file从txt载入矩阵
		~Matrix(); // destruction删除

		void set(int r, int c, T value); // row, column, value
		T get(int r, int c); // row, column
		int rows(); // number of rows
		int columns(); // number of columns

		void print(); // print the matrix
		Matrix * copy(); // copy itself to a new matrix

		void saveTxt(const char * filename); // save matrix to txt file

		// B=M'转置
		Matrix * transpose();
		// B=M(r1:r2,c1:c2)
		Matrix * sub(int r1, int r2, int c1, int c2); // submatrix子矩阵
		// B=|M|
		Matrix * abs(); // absolute values取绝对值

		// numbers of matrix
		T trace(); // trace
		double fnorm(); // Frobenius norm
		double pnorm(double p); // p-norm
		T maxEl(int &r, int &c); // max element
		T minEl(int &r, int &c); // min element
		double mean(); // mean of elements
		T sum(); // sum of elements
		double std(); // standard deviation of elements标准差


		// M=M+a
		void addNumberSelf(T value); // add a number to itself in space
		// M=M*a
		void multiplyNumberSelf(T value); // add a number to itself in space

		// M=M+A
		void addMatrixSelf(Matrix * A); // add a matrix to itself in space
		// M=M.*A
		void dotMultiplyMatrixSelf(Matrix * A); // dot multiply a matrix to itself in space

		// B=M+A
		Matrix * addMatrixNew(Matrix * A); // add a matrix to itself with new matrix
		// B=M.*A
		Matrix * dotMultiplyMatrixNew(Matrix * A); // dot multiply a matrix to itself with new matrix
		// B=M*A
		Matrix * multiplyMatrixNew(Matrix * A); // multiply a matrix to itself with new matrix

		// Multidimensional scaling (MDS)
		// This function re-implements Laurens van der Maaten's MDS in his Matlab Toolbox for Dimensionality Reduction
		// The Matlab MDS can be downloaded at http://crcv.ucf.edu/source/dimension
		Matrix<double> * MDS_UCF(int dim, int iter);

	private:
		int rows_;
		int columns_;
		T **v;
	};



	/**********************************************
	* Utilities part
	**********************************************/

	template<class T>
	T min(T v1, T v2)
	{
		if (v1 < v2) return v1;
		else return v2;
	}

	template<class T>
	T max(T v1, T v2)
	{
		if (v1 > v2) return v1;
		else return v2;
	}

	template<class T>
	void swap(T &v1, T &v2)
	{
		T v3 = v1;
		v1 = v2;
		v2 = v3;
	}

	template<class T>
	double sign(T v)
	{
		if (v > 0) return 1.0;
		else if (v < 0) return -1.0;
		else return 0.0;
	}

	/**********************************************
	* Implementation part
	**********************************************/

	template<class T>
	Matrix<T>::Matrix(int rows, int columns) // initialization without assigning values
	{
		if (rows < 1 || columns < 1)
		{
			printf("Invalid construction arguments: rows=%d, columns=%d\n", rows, columns);
			exit(1);
		}

		rows_ = rows;
		columns_ = columns;

		v = new T *[rows];
		for (int i = 0; i < rows; i++)
		{
			v[i] = new T[columns];
		}
	}

	template<class T>
	Matrix<T>::Matrix(int rows, int columns, T value) // initialization with all same values
	{
		if (rows < 1 || columns < 1)
		{
			printf("Invalid construction arguments: rows=%d, columns=%d\n", rows, columns);
			exit(1);
		}

		rows_ = rows;
		columns_ = columns;

		v = new T *[rows];
		for (int i = 0; i < rows; i++)
		{
			v[i] = new T[columns];

			for (int j = 0; j < columns; j++)
			{
				v[i][j] = value;
			}
		}
	}

	template<class T>
	Matrix<T>::Matrix(int rows, int columns, std::string type) // special matrix such as I
	{
		if (rows < 1 || columns < 1)
		{
			printf("Invalid construction arguments: rows=%d, columns=%d\n", rows, columns);
			exit(1);
		}
		rows_ = rows;
		columns_ = columns;

		v = new T *[rows];
		for (int i = 0; i < rows; i++)
		{
			v[i] = new T[columns];
		}

		if (type.compare("I") == 0)
		{
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < columns; j++)
				{
					if (i == j) v[i][j] = (T)1;
					else v[i][j] = (T)0;
				}
			}
		}

		else if (type.compare("rand") == 0) // all elements between 0 and 1
		{
			srand(time(NULL));
			int r1;
			double r2;
			for (int i = 0; i < rows_; i++)
			{
				for (int j = 0; j < columns_; j++)
				{
					r1 = rand()*rand() + rand()*rand() + rand();
					if (r1 < 0) r1 = -r1;
					r2 = double(r1 % 1000001) / 1000000;

					v[i][j] = (T)r2;
				}
			}
		}

		else if (type.compare("rand_int") == 0)
		{
			srand(time(NULL));
			for (int i = 0; i < rows_; i++)
			{
				for (int j = 0; j < columns_; j++)
				{
					v[i][j] = (T)rand();
				}
			}
		}

		else if (type.compare("randperm") == 0) // random permutation, each column is a randperm vector of size rows*1
		{
			srand(time(NULL));
			for (int j = 0; j < columns; j++)
			{
				for (int i = 0; i < rows; i++)
				{
					v[i][j] = i + 1;
				}

				for (int i = 0; i < rows; i++)
				{
					int k = rand() % rows;
					if (k >= rows || k < 0)
					{
						printf("Invalid row index: %d\n", k);
						exit(1);
					}
					T temp = v[i][j];
					v[i][j] = v[k][j];
					v[k][j] = temp;
				}
			}
		}

		else
		{
			printf("Undefined matrix type: %s\n", type.c_str());
			exit(1);
		}
	}

	template<class T>
	Matrix<T>::Matrix(const char * filename)//read matrix from txt
	{
		FILE * pFile;
		// first pass: matrix size
		int rows = 0;
		int columns = 0;

		pFile = fopen(filename, "r");
		if (pFile == NULL)
		{
			printf("File \"%s\" cannot be found.\n", filename);
			exit(1);
		}
		char line[MAX_LINE_LENGTH];
		char * token = NULL;
		while (fgets(line, MAX_LINE_LENGTH, pFile) != NULL)
		{
			rows++;
			if (rows == 1) // count the number of columns
			{
				token = strtok(line, " ,\t");
				while (token != NULL && token[0] >= 32)
				{
					columns++;
					token = strtok(NULL, " ,\t");
				}
			}
			else // check whether every row contains the same number of elements with the first row
			{
				int check = 0;
				token = strtok(line, " ,\t");
				while (token != NULL && token[0] >= 32)
				{
					check++;
					token = strtok(NULL, " ,\t");
				}
				if (check < columns)
				{
					rows--;
					break;
				}
			}
		}
		fclose(pFile);
		printf("Reading matrix from file \"%s\": %d rows, %d columns\n", filename, rows, columns);

		// second pass: read data
		rows_ = rows;
		columns_ = columns;
		v = new T *[rows];
		for (int i = 0; i < rows; i++)
		{
			v[i] = new T[columns];
		}

		pFile = fopen(filename, "r");
		if (pFile == NULL)
		{
			printf("File \"%s\" cannot be found.\n", filename);
			exit(1);
		}
		int i = 0;
		while (fgets(line, MAX_LINE_LENGTH, pFile) != NULL)
		{
			if (i >= rows) break;
			for (int j = 0; j < columns; j++)
			{
				if (j == 0) token = strtok(line, " ,\t");
				else token = strtok(NULL, " ,\t");
				v[i][j] = (T)atof(token);
			}
			i++;
		}
		fclose(pFile);
	}

	template<class T>
	Matrix<T>::~Matrix() // destruction
	{
		for (int i = 0; i < rows_; i++)
		{
			delete[](T *)v[i];
		}
		delete[] v;
	}

	template<class T>
	void Matrix<T>::set(int r, int c, T value) // set matrix row i, column j, as value
	{
		if (r < 0 || r >= rows_ || c < 0 || c >= columns_)
		{
			printf("Invalid index in set(): r=%d, c=%d\n", r, c);
			exit(1);
		}
		v[r][c] = value;
	}

	template<class T>
	T Matrix<T>::get(int r, int c) //get matrix row i, column j
	{
		if (r < 0 || r >= rows_ || c < 0 || c >= columns_)
		{
			printf("Invalid index in get(): r=%d, c=%d\n", r, c);
			exit(1);
		}
		return v[r][c];
	}

	template<class T>
	int Matrix<T>::rows() // number of rows
	{
		return rows_;
	}

	template<class T>
	int Matrix<T>::columns() // number of columns
	{
		return columns_;
	}

	template<class T>
	void Matrix<T>::print() // print the matrix
	{
		printf("\n");
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				printf("%8.5f    ", (double)v[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}

	template<class T>
	Matrix<T> * Matrix<T>::copy() // copy itself to a new matrix
	{
		Matrix<T> * A = new Matrix<T>(rows_, columns_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				A->set(i, j, v[i][j]);
			}
		}
		return A;
	}

	template<class T>
	void Matrix<T>::saveTxt(const char * filename)//save matrix as file
	{
		FILE * pFile;
		pFile = fopen(filename, "w");
		if (pFile == NULL)
		{
			printf("Cannot save to file \"%s\".\n", filename);
			exit(1);
		}
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				fprintf(pFile, "%f  ", (double)v[i][j]);
			}
			fprintf(pFile, "\n");
		}
		fclose(pFile);
		printf("Matrix saved to file \"%s\"\n", filename);
	}

	template<class T>
	Matrix<T> * Matrix<T>::transpose()
	{
		Matrix<T> * A = new Matrix<T>(columns_, rows_);
		for (int i = 0; i < columns_; i++)
		{
			for (int j = 0; j < rows_; j++)
			{
				A->set(i, j, v[j][i]);
			}
		}
		return A;
	}

	template<class T>
	Matrix<T> * Matrix<T>::sub(int r1, int r2, int c1, int c2) // submatrix
	{
		if (r1 < 0 || r1 >= rows_ || r2 < 0 || r2 >= rows_ || r2 < r1 || c1 < 0 || c1 >= columns_ || c2<0 || c2>columns_ || c2 < c1)
		{
			printf("Invalid submatrix indices.\n");
			exit(1);
		}

		int newRows = r2 - r1 + 1;
		int newColumns = c2 - c1 + 1;
		Matrix<T> * A = new Matrix<T>(newRows, newColumns);
		for (int i = 0; i < newRows; i++)
		{
			for (int j = 0; j < newColumns; j++)
			{
				A->set(i, j, v[i + r1][j + c1]);
			}
		}
		return A;
	}

	template<class T>
	Matrix<T> * Matrix<T>::abs() // absolute values
	{
		Matrix<T> * A = new Matrix<T>(rows_, columns_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				A->set(i, j, v[i][j] > 0 ? v[i][j] : -v[i][j]);
			}
		}
		return A;
	}

	template<class T>
	T Matrix<T>::trace() // trace   对角值
	{
		T x = 0;
		for (int i = 0; i < min<int>(rows_, columns_); i++)
		{
			x += v[i][i];
		}
		return x;
	}

	template<class T>
	double Matrix<T>::fnorm() // Frobenius norm
	{
		double x = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				x += (v[i][j] * v[i][j]);
			}
		}
		return sqrt(x);
	}

	template<class T>
	double Matrix<T>::pnorm(double p) // p-norm
	{
		double x = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				x += pow(fabs((double)v[i][j]), p);
			}
		}
		return pow(x, 1 / p);
	}

	template<class T>
	T Matrix<T>::maxEl(int &r, int &c) // max element
	{
		T x = v[0][0];
		r = 0;
		c = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				if (v[i][j] > x)
				{
					x = v[i][j];
					r = i;
					c = j;
				}
			}
		}
		return x;
	}

	template<class T>
	T Matrix<T>::minEl(int &r, int &c) // min element
	{
		T x = v[0][0];
		r = 0;
		c = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				if (v[i][j] < x)
				{
					x = v[i][j];
					r = i;
					c = j;
				}
			}
		}
		return x;
	}

	template<class T>
	double Matrix<T>::mean() // mean of elements
	{
		double x = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				x += (double)v[i][j];
			}
		}
		return x / rows_ / columns_;
	}

	template<class T>
	T Matrix<T>::sum() // sum of elements
	{
		T x = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				x += v[i][j];
			}
		}
		return x;
	}

	template<class T>
	double Matrix<T>::std() // standard deviation of elements
	{
		double m = mean();
		double s = 0;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				s += (v[i][j] - m)*(v[i][j] - m);
			}
		}
		s = (s / rows_) / columns_;
		return sqrt(s);
	}

	template<class T>
	void Matrix<T>::addNumberSelf(T value) // add a number to itself in space
	{
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				v[i][j] += value;
			}
		}
	}

	template<class T>
	void Matrix<T>::multiplyNumberSelf(T value) // multi a number to itself in space
	{
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				v[i][j] *= value;
			}
		}
	}

	template<class T>
	void Matrix<T>::addMatrixSelf(Matrix * A) // add a matrix to itself in space
	{
		if (rows_ != A->rows() || columns_ != A->columns())
		{
			printf("Unmatched matrix sizes in matrix summation.\n");
			exit(1);
		}

		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				v[i][j] += A->get(i, j);
			}
		}
	}

	template<class T>
	void Matrix<T>::dotMultiplyMatrixSelf(Matrix * A) // dot multiply a matrix to itself in space
	{
		if (rows_ != A->rows() || columns_ != A->columns())
		{
			printf("Unmatched matrix sizes in matrix dot multiplication.\n");
			exit(1);
		}

		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				v[i][j] *= A->get(i, j);
			}
		}
	}

	template<class T>
	Matrix<T> * Matrix<T>::addMatrixNew(Matrix * A) // add a matrix to itself with new matrix
	{
		if (rows_ != A->rows() || columns_ != A->columns())
		{
			printf("Unmatched matrix sizes in matrix summation.\n");
			exit(1);
		}

		Matrix<T> * B = new Matrix<T>(rows_, columns_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				B->set(i, j, v[i][j] + A->get(i, j));
			}
		}
		return B;
	}

	template<class T>
	Matrix<T> * Matrix<T>::dotMultiplyMatrixNew(Matrix * A) // dot multiply a matrix to itself with new matrix
	{
		if (rows_ != A->rows() || columns_ != A->columns())
		{
			printf("Unmatched matrix sizes in matrix dot multiplication.\n");
			exit(1);
		}

		Matrix<T> * B = new Matrix<T>(rows_, columns_);
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < columns_; j++)
			{
				B->set(i, j, v[i][j] * A->get(i, j));
			}
		}
		return B;
	}

	template<class T>
	Matrix<T> * Matrix<T>::multiplyMatrixNew(Matrix * A) // multiply a matrix to itself with new matrix
	{
		if (columns_ != A->rows())
		{
			printf("Unmatched matrix sizes in matrix multiplication.\n");
			exit(1);
		}

		Matrix<T> * B = new Matrix<T>(rows_, A->columns());
		T temp;
		for (int i = 0; i < rows_; i++)
		{
			for (int j = 0; j < A->columns(); j++)
			{
				temp = 0;
				for (int k = 0; k < columns_; k++)
				{
					temp += (v[i][k] * A->get(k, j));
				}
				B->set(i, j, temp);
			}
		}
		return B;
	}

	template<class T>
	Matrix<double> * Matrix<T>::MDS_UCF(int dim, int iter)
	{
		if (rows() != columns()) { printf("Input distance matrix to MDS is not square.\n"); exit(1); }
		if (dim < 1) { printf("Invalid dimension for MDS.\n"); exit(1); }
		if (iter < 1) { printf("Invalid number of iterations for MDS.\n"); exit(1); }

		Matrix<double> * X = new Matrix<double>(rows(), dim, "rand");
		double D_mean = mean(); // mean value of distance matrix
		X->addNumberSelf(-0.5); // move to the center
		X->multiplyNumberSelf(0.1*D_mean / (1.0 / 3.0*sqrt((double)dim))); // before this step, mean distance is 1/3*sqrt(d)

		double lr = 0.05; // learning rate
		double r = 2; // metric
		int n = rows(); // number of vectors

		Matrix<double> * dh = new Matrix<double>(n, n, 0.0);
		Matrix<double> * pmat = new Matrix<double>(n, dim);
		Matrix<double> * dhdum = new Matrix<double>(n, 1);
		Matrix<double> * dhmat = new Matrix<double>(n - 1, dim, 0);

		Matrix<int> * RP = new Matrix<int>(n, iter, "randperm"); // the matrix for random permutation numbers
		int i, j;
		double temp;
		int m;

		//printf("MDS iteration:");
		for (int it = 0; it < iter; it++) // iterations
		{
			//if (it % 10 == 0) printf("\n");
			//printf("%3d  ", it + 1);
			for (int rp = 0; rp < n; rp++) // work on each vector in a randomly permuted order
			{
				m = RP->get(rp, it) - 1;

				for (i = 0; i < n; i++)
				{
					for (j = 0; j < dim; j++)
					{
						pmat->set(i, j, X->get(m, j) - X->get(i, j));
					}
				}

				for (i = 0; i < n; i++)
				{
					temp = 0;
					for (j = 0; j < dim; j++)
					{
						temp += pow(fabs(pmat->get(i, j)), r);
					}
					dhdum->set(i, 0, pow(temp, 1 / r));
				}

				for (i = 0; i < n; i++)
				{
					if (i == m) continue;

					dh->set(m, i, dhdum->get(i, 0));
					dh->set(i, m, dhdum->get(i, 0));
				}

				for (i = 0; i < n - 1; i++)
				{
					int ii = i;
					if (i >= m) ii = i + 1;
					temp = lr * (dhdum->get(ii, 0) - get(ii, m)) * pow(dhdum->get(ii, 0), 1 - r);
					for (j = 0; j < dim; j++)
					{
						dhmat->set(i, j, temp);
					}
				}

				for (i = 0; i < n - 1; i++)
				{
					int ii = i;
					if (i >= m) ii = i + 1;
					for (j = 0; j < dim; j++)
					{
						temp = X->get(ii, j);
						temp += dhmat->get(i, j) * pow(fabs(pmat->get(ii, j)), r - 1) * sign<double>(pmat->get(ii, j));

						X->set(ii, j, temp);
					}
				}
			}
		}

		delete dh;
		delete pmat;
		delete dhdum;
		delete dhmat;
		delete RP;

		return X;
	}

	/**********************************************
	* Algorithm part
	**********************************************/

	// Calculate the pairwise interpoint Euclidean distances
	// X is data matrix, D is distance matrix
	/*void EuclideanDistanceMatrix(Matrix<double> * X, Matrix<double> * D)
	{
	int i, j, k;
	double temp;
	if (D == NULL)
	{
	printf("Input matrix pointer is NULL.\n");
	exit(1);
	}
	else if (X->rows() != D->rows() || X->rows() != D->columns())
	{
	printf("Invalid distance matrix dimension.\n");
	exit(1);
	}
	for (i = 0; i<D->rows(); i++) D->set(i, i, 0.0);
	for (i = 0; i<D->rows() - 1; i++)
	{
	for (j = i + 1; j<D->columns(); j++)
	{
	temp = 0;
	for (k = 0; k<X->columns(); k++)
	{
	temp += pow(X->get(i, k) - X->get(j, k), 2);
	}
	D->set(i, j, sqrt(temp));
	}
	}
	for (i = 1; i<D->rows(); i++)
	{
	for (j = 0; j<i; j++)
	{
	D->set(i, j, D->get(j, i));
	}
	}
	}*/

	// Copy all elements of X to Y
	/*void ElementCopy(Matrix<double> * X, Matrix<double> * Y)
	{
	if (Y == NULL)
	{
	printf("Input matrix pointer is NULL.\n");
	exit(1);
	}
	else if (X->rows() != Y->rows() || X->columns() != Y->columns())
	{
	printf("Invalid matrix dimension.\n");
	exit(1);
	}
	for (int i = 0; i<X->rows(); i++)
	{
	for (int j = 0; j<X->columns(); j++)
	{
	Y->set(i, j, X->get(i, j));
	}
	}
	}*/

	// Multidimensional scaling (MDS)
	// This function re-implements Laurens van der Maaten's MDS in his Matlab Toolbox for Dimensionality Reduction
	// The Matlab MDS can be downloaded at http://crcv.ucf.edu/source/dimension
	//Matrix<double> * MDS_UCF(Matrix<double> * D, Matrix<double> * X0, int dim, int iter)
	//{
	//	if (D->rows() != D->columns())
	//	{
	//		printf("Input distance matrix to MDS is not square.\n");
	//		exit(1);
	//	}
	//	if (dim<1)
	//	{
	//		printf("Invalid dimension for MDS.\n");
	//		exit(1);
	//	}
	//	if (iter<1)
	//	{
	//		printf("Invalid number of iterations for MDS.\n");
	//		exit(1);
	//	}

	//	Matrix<double> * X = NULL;

	//	// with initialization
	//	if (X0 != NULL)
	//	{
	//		if (X0->rows() != D->rows() || X0->columns() != dim)
	//		{
	//			printf("Input initialization to MDS has invalid dimension.\n");
	//			exit(1);
	//		}
	//		X = X0->copy();
	//	}
	//	// without initialization
	//	else
	//	{
	//		X = new Matrix<double>(D->rows(), dim, "rand");
	//		double D_mean = D->mean(); // mean value of distance matrix
	//		X->addNumberSelf(-0.5); // move to the center
	//		X->multiplyNumberSelf(0.1*D_mean / (1.0 / 3.0*sqrt((double)dim))); // before this step, mean distance is 1/3*sqrt(d)
	//	}

	//	double lr = 0.05; // learning rate
	//	double r = 2; // metric
	//	int n = D->rows(); // number of vectors


	//	Matrix<double> * dh = new Matrix<double>(n, n, 0.0);
	//	Matrix<double> * pmat = new Matrix<double>(n, dim);
	//	Matrix<double> * dhdum = new Matrix<double>(n, 1);
	//	Matrix<double> * dhmat = new Matrix<double>(n - 1, dim, 0);

	//	Matrix<int> * RP = new Matrix<int>(n, iter, "randperm"); // the matrix for random permutation numbers
	//	int i, j;
	//	double temp;
	//	int m;

	//	printf("MDS iteration:");
	//	for (int it = 0; it<iter; it++) // iterations
	//	{
	//		if (it % 10 == 0) printf("\n");
	//		printf("%3d  ", it + 1);
	//		for (int rp = 0; rp<n; rp++) // work on each vector in a randomly permuted order
	//		{
	//			m = RP->get(rp, it) - 1;

	//			for (i = 0; i<n; i++)
	//			{
	//				for (j = 0; j<dim; j++)
	//				{
	//					pmat->set(i, j, X->get(m, j) - X->get(i, j));
	//				}
	//			}

	//			for (i = 0; i<n; i++)
	//			{
	//				temp = 0;
	//				for (j = 0; j<dim; j++)
	//				{
	//					temp += pow(fabs(pmat->get(i, j)), r);
	//				}
	//				dhdum->set(i, 0, pow(temp, 1 / r));
	//			}

	//			for (i = 0; i<n; i++)
	//			{
	//				if (i == m) continue;

	//				dh->set(m, i, dhdum->get(i, 0));
	//				dh->set(i, m, dhdum->get(i, 0));
	//			}

	//			for (i = 0; i<n - 1; i++)
	//			{
	//				int ii = i;
	//				if (i >= m) ii = i + 1;
	//				temp = lr * (dhdum->get(ii, 0) - D->get(ii, m)) * pow(dhdum->get(ii, 0), 1 - r);
	//				for (j = 0; j<dim; j++)
	//				{
	//					dhmat->set(i, j, temp);
	//				}
	//			}

	//			for (i = 0; i<n - 1; i++)
	//			{
	//				int ii = i;
	//				if (i >= m) ii = i + 1;
	//				for (j = 0; j<dim; j++)
	//				{
	//					temp = X->get(ii, j);
	//					temp += dhmat->get(i, j) * pow(fabs(pmat->get(ii, j)), r - 1) * sign<double>(pmat->get(ii, j));

	//					X->set(ii, j, temp);
	//				}
	//			}
	//		}
	//	}

	//	printf("\n");

	//	delete dh;
	//	delete pmat;
	//	delete dhdum;
	//	delete dhmat;
	//	delete RP;

	//	return X;
	//}

	// Multidimensional scaling (MDS) with SMACOF
	// This code re-implements Michael Bronstein's SMACOF in his Matlab Toolbox for Surface Comparison and Analysis
	// The Matlab SMACOF can be downloaded at http://tosca.cs.technion.ac.il/
	//Matrix<double> * MDS_SMACOF(Matrix<double> * D, Matrix<double> * X0, int dim, int iter)
	//{
	//	if (D->rows() != D->columns())
	//	{
	//		printf("Input distance matrix to MDS is not square.\n");
	//		exit(1);
	//	}
	//	if (dim<1)
	//	{
	//		printf("Invalid dimension for MDS.\n");
	//		exit(1);
	//	}
	//	if (iter<1)
	//	{
	//		printf("Invalid number of iterations for MDS.\n");
	//		exit(1);
	//	}

	//	Matrix<double> * X = NULL;

	//	// with initialization
	//	if (X0 != NULL)
	//	{
	//		if (X0->rows() != D->rows() || X0->columns() != dim)
	//		{
	//			printf("Input initialization to MDS has invalid dimension.\n");
	//			exit(1);
	//		}
	//		X = X0->copy();
	//	}
	//	// without initialization
	//	else
	//	{
	//		X = new Matrix<double>(D->rows(), dim, "rand");
	//		double D_mean = D->mean(); // mean value of distance matrix
	//		X->addNumberSelf(-0.5); // move to the center
	//		X->multiplyNumberSelf(0.1*D_mean / (1.0 / 3.0*sqrt((double)dim))); // before this step, mean distance is 1/3*sqrt(d)
	//	}


	//	Matrix<double> * Z = X->copy();
	//	Matrix<double> * D_ = new Matrix<double>(D->rows(), D->columns());
	//	Matrix<double> * B = new Matrix<double>(D->rows(), D->columns());
	//	int i, j, k;
	//	double temp;

	//	EuclideanDistanceMatrix(X, D_);

	//	printf("MDS iteration:");
	//	for (int it = 0; it<iter; it++) // iterations
	//	{
	//		if (it % 10 == 0) printf("\n");
	//		printf("%3d  ", it + 1);

	//		// B = calc_B(D_,D);
	//		for (i = 0; i<D->rows(); i++)
	//		{
	//			for (j = 0; j<D->columns(); j++)
	//			{
	//				if (i == j || fabs(D_->get(i, j))<EPSILON)
	//				{
	//					B->set(i, j, 0.0);
	//				}
	//				else
	//				{
	//					B->set(i, j, -D->get(i, j) / D_->get(i, j));
	//				}
	//			}
	//		}

	//		for (j = 0; j<D->columns(); j++)
	//		{
	//			temp = 0;
	//			for (i = 0; i<D->rows(); i++)
	//			{
	//				temp += B->get(i, j);
	//			}
	//			B->set(j, j, -temp);
	//		}

	//		// X = B*Z/size(D,1);
	//		for (i = 0; i<X->rows(); i++)
	//		{
	//			for (j = 0; j<X->columns(); j++)
	//			{
	//				temp = 0;
	//				for (k = 0; k<B->columns(); k++)
	//				{
	//					temp += (B->get(i, k)*Z->get(k, j));
	//				}
	//				X->set(i, j, temp / (double)D->rows());
	//			}
	//		}

	//		// D_ = calc_D (X);
	//		EuclideanDistanceMatrix(X, D_);

	//		// Z = X;
	//		ElementCopy(X, Z);
	//	}

	//	printf("\n");

	//	delete Z;
	//	delete D_;
	//	delete B;

	//	return X;
	//}
}

#endif

//------------------------



struct node {
	string feature;
	vector<int> ve;
}tree[630];
map<string, int> s2i;
map<int, string> i2s;
int cnt = 0;
int num_feature = 0;
vector<string> features;
void build(int root, const vector<string> vec, int k) {
	if (k == vec.size()) {
		tree[root].feature = features[num_feature];
		return;
	}
	if (s2i.find(vec[k]) == s2i.end()) {
		s2i[vec[k]] = ++cnt;
		i2s[cnt] = vec[k];
		tree[root].ve.push_back(cnt);
		build(cnt, vec, k + 1);
	}
	else {
		build(s2i[vec[k]], vec, k + 1);
	}
}

void print(int rt) {
	cout << i2s[rt] << " ";
	for (int i = 0;i < tree[rt].ve.size();i++) {
		print(tree[rt].ve[i]);
	}
	cout << endl;
}

std::string trim(std::string currStr)
{
	if (currStr.empty()) return currStr;

	currStr.erase(0, currStr.find_first_not_of(" "));
	currStr.erase(currStr.find_last_not_of(" ") + 1);
	return currStr;
}

void split(const std::string& currStr, std::string delim, std::vector<std::string > & ret)
{
	size_t last = 0;
	size_t index = currStr.find_first_of(delim, last);

	while (index != std::string::npos)
	{
		ret.push_back(trim(currStr.substr(last, index - last)));
		last = index + 1;
		index = currStr.find_first_of(delim, last);
	}
	if (index - last > 0)  ret.push_back(trim(currStr.substr(last, index - last)));
}

std::string replace_all(std::string currStr, std::string a, std::string b)
{
	if (currStr.empty()) return currStr;

	int len = currStr.length();
	for (int pos = 0; pos < len; ++pos)
		if (currStr.substr(pos, 1) == a)
			currStr.replace(pos, 1, b);
	return currStr;
}


std::string getFileName(const std::string & str_arg_path)
{
	std::string base = str_arg_path.substr(str_arg_path.find_last_of("/\\") + 1);
	const size_t period_idx = base.rfind('.');
	if (std::string::npos != period_idx) base.erase(period_idx);
	return base;
}

bool endsWith(std::string const & currStr, std::string const & ending)
{
	if (ending.size() > currStr.size()) return false;
	return std::equal(ending.rbegin(), ending.rend(), currStr.rbegin());
}

void nextnode(boost::property_tree::ptree currnode, std::string child, std::vector<std::string > & ret)
{
	boost::property_tree::ptree temp;
	temp = currnode.get_child(child);
	for (boost::property_tree::ptree::iterator it = currnode.begin(); it != currnode.end(); ++it)
	{
		string key = it->first;//key ID
		ret.push_back(trim(key));
	}
}

int main(int argc, char* argv[])
{
	std::vector<std::string> vec_embeddingFiles, vec_data, vec_top, vec_embed, curnum, vec_names, vec_info, vec_info2, tmp_name, tmp_node, vec_db, vec_dis, in_embed, vec_10, vec_seq, vec_noseq;
	std::string str_dbName, str_DataEmbedFile, str_outPutFile, child1, curdbase, str_dbDir, str_queryEmbedFile, turelab, alldiff, treename, jjj, outinfoall;

	smat::Matrix<double> * Mat_dist = new smat::Matrix<double>(11, 11, "rand");
	smat::Matrix<double> * Mat_pcoa = new smat::Matrix<double>(11, 2, "rand");


	//Mat_dist->set(0, 0, 0);
		//Matrix<double> * MDS_UCF(2, int iter);

		//float wordVec[256][10],contextVec[256][10],wordBias[256],contextBias[256];
	int i_k = 4, i_dim = 20, i_genomeCnt = 0;
	boost::property_tree::ptree datum;
	clock_t startTime, endTime;
	startTime = clock();
	bool readEmbedSuccess = false;




	for (int i = 1; i < argc; ++i)
	{
		if (!strcmp(argv[i], "-I") || !strcmp(argv[i], "-i")) str_queryEmbedFile = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-K") || !strcmp(argv[i], "-k")) i_k = atoi(argv[++i]);
		//else if (!strcmp(argv[i], "-M") || !strcmp(argv[i], "-m")) treename = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-D") || !strcmp(argv[i], "-d")) str_dbDir = std::string(argv[++i]);
		else if (!strcmp(argv[i], "-S") || !strcmp(argv[i], "-s")) str_dbName = std::string(argv[++i]);
		//else if (!strcmp(argv[i], "-T") || !strcmp(argv[i], "-t")) i_top = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-DIM") || !strcmp(argv[i], "-dim")) i_dim = atoi(argv[++i]);
		else if (!strcmp(argv[i], "-M") || !strcmp(argv[i], "-m")) treename = std::string(argv[++i]);
	}

	printf("Start initializing the parameters... \n");
	unsigned long vocabLen = (unsigned long)pow(BASE, i_k);
	REAL* wordBias = new REAL[vocabLen]; memset(wordBias, 0, sizeof(REAL) * vocabLen); //创建长为256的数组，并给所有元素赋初值为0
	REAL* contextBias = new REAL[vocabLen]; memset(contextBias, 0, sizeof(REAL) * vocabLen);
	REAL** wordVec = new REAL*[vocabLen]; REAL** contextVec = new REAL*[vocabLen];
	REAL** srcProbTransMat = new REAL*[vocabLen]; REAL** trgtProbTransMat = new REAL*[vocabLen];
	split(str_dbName, ",", vec_db);
	for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx)
	{
		wordVec[vocabIdx] = new REAL[i_dim]; memset(wordVec[vocabIdx], 0, sizeof(REAL) * i_dim);
		contextVec[vocabIdx] = new REAL[i_dim]; memset(contextVec[vocabIdx], 0, sizeof(REAL) * i_dim);

		srcProbTransMat[vocabIdx] = new REAL[vocabLen]; memset(srcProbTransMat[vocabIdx], 0, sizeof(REAL) * vocabLen);
		trgtProbTransMat[vocabIdx] = new REAL[vocabLen]; memset(trgtProbTransMat[vocabIdx], 0, sizeof(REAL) * vocabLen);
	}
	REAL* probVec = new REAL[vocabLen]; memset(probVec, 0, sizeof(REAL) * vocabLen);



	printf("Start Loading the embedding file... \n");
	std::ifstream embedStream(str_queryEmbedFile.c_str());
	if (embedStream.is_open())
	{
		std::string str_tmp_line = ""; unsigned long currVocabIdx = 0;
		while (getline(embedStream, str_tmp_line))
		{
			if (str_tmp_line.empty() || "" == str_tmp_line) continue;
			if (str_tmp_line.substr(0, 1) == ">")
			{
				currVocabIdx = 0;
				memset(wordBias, 0, sizeof(REAL) * vocabLen); memset(contextBias, 0, sizeof(REAL) * vocabLen);
				for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx)
				{
					memset(wordVec[vocabIdx], 0, sizeof(REAL) * i_dim);
					memset(contextVec[vocabIdx], 0, sizeof(REAL) * i_dim);
					memset(srcProbTransMat[vocabIdx], 0, sizeof(REAL) * vocabLen);
				}
				turelab = str_tmp_line;

				continue;
			}
			vec_embed.clear();
			split(str_tmp_line, ",", vec_embed);

			for (int j = 0; j < i_dim; ++j) { wordVec[currVocabIdx][j] = atof(vec_embed[j].c_str()); contextVec[currVocabIdx][j] = atof(vec_embed[j + i_dim + 1].c_str()); }
			wordBias[currVocabIdx] = atof(vec_embed[i_dim].c_str()); contextBias[currVocabIdx] = atof(vec_embed[2 * i_dim + 1].c_str());
			currVocabIdx++;

			if (currVocabIdx == vocabLen)
			{
				for (unsigned long row = 0; row < vocabLen; ++row)
				{
					memset(probVec, 0, sizeof(REAL) * vocabLen); REAL colSum = 0;
					for (unsigned long col = 0; col < vocabLen; ++col)
					{
						REAL sum = 0;
						for (int i = 0; i < i_dim; ++i) sum += wordVec[row][i] * contextVec[col][i];
						probVec[col] = exp(sum + contextBias[col] + wordBias[row]); colSum += probVec[col];
						//cout<<probVec[col]<<endl;
					}
					for (unsigned long col = 0; col < vocabLen; ++col) srcProbTransMat[row][col] = (probVec[col] / colSum);
				}
				readEmbedSuccess = true; break;
			}

		}
		embedStream.close();
	}
	else
		fprintf(stderr, "[ERROR!] The query embedding file %s doesn't exist!\n", str_queryEmbedFile.c_str());
	if (!readEmbedSuccess) fprintf(stderr, "[ERROR!] The query embedding file is incomplete !\n");

	std::priority_queue<std::pair<REAL, std::string>> distIDQueue;
	boost::property_tree::ptree root;
	boost::property_tree::ptree items, zhong;
	boost::property_tree::ptree item, item1, item2, item3, item4, item5, item6, item7;
	//boost::property_tree::read_json<boost::property_tree::ptree>(str_dbDir, root);
	boost::property_tree::read_json(treename, root);
	items = root;
	int inw = 1;
	int name = 0;
	std::string str_currDbDir = str_dbDir, allID = "", allkey = "";

	std::string str_dbData = str_currDbDir; str_dbData.append(".data"); std::cout << str_dbData << std::endl;
	std::string str_dbInfo = str_currDbDir; str_dbInfo.append(".info1"); std::cout << str_dbInfo << std::endl;

	//std::ifstream ifsData("/home/wangying/wangkun/CRAFT/sampling/bacteria_k4_t1_v10.data", std::ios::in | std::ios::binary);
	std::ifstream ifsData(str_dbData, std::ios::in | std::ios::binary);
	std::ifstream ifsinfo1(str_dbInfo, std::ios::in);
	//std::ifstream ifsData(str_dbData.c_str(), std::ios::in | std::ios::binary);
	std::ifstream ifsInfo(str_dbInfo.c_str()); std::string str_seqName = "";
	int dn = 0;
	child1 = "all";
	if (str_dbName == "All")str_dbName = "";
	vec_top.push_back(":" + str_dbName);
	if (str_dbName != "")dn = 1;
	split(child1, ",", vec_data);


    name = 0;
    for (int c = 0;c < vec_top.size();c++)
    {
        vec_data.clear();
        split(vec_top[c], ":", vec_data);
        item = items;
        for (boost::property_tree::ptree::iterator it = item.begin(); it != item.end(); ++it)
        {

            string key = it->first;//key ID
            item1 = item.get_child(key);//jie
            for (boost::property_tree::ptree::iterator it1 = item1.begin(); it1 != item1.end(); ++it1)
            {
                string key1 = it1->first;

                item2 = item1.get_child(key1);

                for (boost::property_tree::ptree::iterator it2 = item2.begin(); it2 != item2.end(); ++it2)
                {

                    string key2 = it2->first;//key ID


                    item3 = item2.get_child(key2);
                    for (boost::property_tree::ptree::iterator it3 = item3.begin(); it3 != item3.end(); ++it3)
                    {

                        string key3 = it3->first;//key ID

                        item4 = item3.get_child(key3);



                        for (boost::property_tree::ptree::iterator it4 = item4.begin(); it4 != item4.end(); ++it4)
                        {

                            string key4 = it4->first;//key ID

                            item5 = item4.get_child(key4);


                            for (boost::property_tree::ptree::iterator it5 = item5.begin(); it5 != item5.end(); ++it5)
                            {

                                string key5 = it5->first;//key ID

                                item6 = item5.get_child(key5);



                                for (boost::property_tree::ptree::iterator it6 = item6.begin(); it6 != item6.end(); ++it6)
                                {

                                    string key6 = it6->first;//key ID
                                     //if(i<1)
                                    item7 = item6.get_child(key6);
                                    //else {item7 = item6.get_child(vec_data[7]);
                                   // }



                                    for (boost::property_tree::ptree::iterator it7 = item7.begin(); it7 != item7.end(); ++it7)
                                    {

                                        string key7 = it7->first;//key ID
                                        string ID = it7->second.get_value<string>(key6);
                                        allID += "," + ID;
                                        name++;

                                        int douhaonum = 0;
                                        basic_string <char>::size_type index = 0;
                                        while ((index = ID.find(',', index)) != string::npos)
                                        {
                                            index++;
                                            douhaonum++;
                                        }

                                        for (int i = 0;i < douhaonum + 1;i++)allkey += "," + key + ":" + key1 + ":" + key2 + ":" + key3 + ":" + key4 + ":" + key5 + ":" + key6 + ":" + key7;

                                    }
                                }
                            }
                        }
                    }
                }
            }if ( dn == 1)break;
        }
    }

    allID.erase(0, 1);
    vec_embed.clear();
    split(allID, ",", vec_embed);

	//std::cout<<'\n'<<vec_embed.size()<<'\n';
    vec_names.clear();
    split(allkey, ",", vec_names);
    allID = "";
    allkey = "";



		int flag = 0;
		name = 0;
		for (int i = 0;i < vec_embed.size();i++)
		{

			stringstream ss;
			ss << vec_embed[i];
			int num;
			ss >> num;

			ifsData.seekg(sizeof(REAL) * vocabLen * (2 * i_dim + 2) * num, std::ios::beg);

			memset(wordBias, 0, sizeof(REAL) * vocabLen); memset(contextBias, 0, sizeof(REAL) * vocabLen);
			for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx)
			{
				memset(wordVec[vocabIdx], 0, sizeof(REAL) * i_dim);
				memset(contextVec[vocabIdx], 0, sizeof(REAL) * i_dim);
			}
			for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx) ifsData.read((char*)wordVec[vocabIdx], sizeof(REAL)*i_dim);
			for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx) ifsData.read((char*)contextVec[vocabIdx], sizeof(REAL)*i_dim);
			ifsData.read((char*)wordBias, sizeof(REAL)*vocabLen);
			string mp;
			mp = std::to_string(wordBias[0]);
			ifsData.read((char*)contextBias, sizeof(REAL)*vocabLen);

			for (unsigned long row = 0; row < vocabLen; ++row)
			{
				memset(probVec, 0, sizeof(REAL) * vocabLen); REAL colSum = 0;
				for (unsigned long col = 0; col < vocabLen; ++col)
				{
					REAL sum = 0;
					for (int i = 0; i < i_dim; ++i) sum += wordVec[row][i] * contextVec[col][i];
					probVec[col] = exp(sum + contextBias[col] + wordBias[row]); colSum += probVec[col];
				}
				for (unsigned long col = 0; col < vocabLen; ++col) trgtProbTransMat[row][col] = (probVec[col] / colSum);

			}

			REAL diff = 0.0;
			for (unsigned long row = 0; row < vocabLen; ++row)
			{
				for (unsigned long col = 0; col < vocabLen; ++col)
					diff += fabs(trgtProbTransMat[row][col] - srcProbTransMat[row][col]);

			}

			getline(ifsInfo, str_seqName);
			name++;
			vec_names[name] = vec_names[name] + ":" + vec_embed[i];
			distIDQueue.push(std::pair<REAL, std::string>(-diff, vec_names[name]));

		}


		std::cout << "distIDQueue size: " << distIDQueue.size() << std::endl;
		std::string outFileName = "";
		curdbase = replace_all(str_dbDir, "/", "-");
		curdbase = replace_all(curdbase, "\\", "-");
		outFileName.append(str_queryEmbedFile).append("_").append(curdbase).append("_").append(str_dbName).append("_dist.txt");
		std::ofstream ofsDiff(outFileName.c_str(), std::ofstream::out);
		std::cout << outFileName << std::endl;
		if (distIDQueue.empty())
			std::cout << "empty!" << std::endl;
		int count = 0;
		int noclap = 0;
		std::vector<std::string>top10;
		child1 = "";

		vec_top.clear();
		vec_dis.clear();
		alldiff = "";
		split(str_queryEmbedFile, "_", in_embed);
		alldiff.append(in_embed[0]).append(",0");

		while (!distIDQueue.empty())
		{
			REAL diff = distIDQueue.top().first; std::string seq = distIDQueue.top().second;
			ofsDiff << "diff:" << diff << "\tid:" << seq << std::endl;
			if (count == 0)
			{
				child1.append(":").append(seq);
			}
			seq = ":" + seq;





            if (count < 10) {
                split(seq, ":", vec_data);
                curnum.push_back(vec_data[vec_data.size() - 1].c_str());
                stringstream ss1;
                ss1 << -diff;
                string currdiff = "";
                ss1 >> currdiff;
                vec_dis.push_back(currdiff);
                Mat_dist->set(0, count + 1, -diff);
                Mat_dist->set(count + 1, 0, -diff);

            }

            Mat_dist->set(0, 0, 0);


            if (count < 1) {
                split(vec_top[0], ":", vec_data);
                stringstream ss;
                ss << vec_data[vec_data.size() - 1];
                int frist;
                ss >> frist;
                ifsinfo1.seekg(0, std::ios::beg);
                for (int k = 0;k < frist + 1;k++)
                    getline(ifsinfo1, str_seqName);
            }

			distIDQueue.pop();
			count++;
		}
		ofsDiff.close();


	std::string outCntName = "";
	outCntName.append(str_queryEmbedFile).append("_").append(curdbase).append("_").append(str_dbName).append(".infoAll");
	std::ofstream outinfoFile(outCntName.c_str(), std::ios::out);
	outinfoall = "dist,pcoaX,pcoaY,id,phylum,class1,order,family,genus,species,infraspecific_name";
	outinfoFile << outinfoall << endl;
	outCntName = "";
	outCntName.append(str_queryEmbedFile).append("_").append(curdbase).append("_").append(str_dbName).append(".infoMat");
	std::ofstream outMatFile(outCntName.c_str(), std::ios::out);

	for (int M = 0;M < vec_dis.size();M++)
	{
		alldiff.append(",").append(vec_dis[vec_dis.size() - 1 - M]);
	}

	outMatFile << alldiff << endl;
	//1to1  top10
	flag = 0;
	name = 0;
	for (int i = 0;i < curnum.size();i++)
	{
		stringstream ss;
		ss << curnum[curnum.size() - 1 - i];
		int num1;
		ss >> num1;
		ifsData.seekg(sizeof(REAL) * vocabLen * (2 * i_dim + 2) * num1, std::ios::beg);

		memset(wordBias, 0, sizeof(REAL) * vocabLen); memset(contextBias, 0, sizeof(REAL) * vocabLen);
		for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx)
		{
			memset(wordVec[vocabIdx], 0, sizeof(REAL) * i_dim);
			memset(contextVec[vocabIdx], 0, sizeof(REAL) * i_dim);
		}
		for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx) ifsData.read((char*)wordVec[vocabIdx], sizeof(REAL)*i_dim);
		for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx) ifsData.read((char*)contextVec[vocabIdx], sizeof(REAL)*i_dim);
		ifsData.read((char*)wordBias, sizeof(REAL)*vocabLen);
		ifsData.read((char*)contextBias, sizeof(REAL)*vocabLen);

		for (unsigned long row = 0; row < vocabLen; ++row)
		{
			memset(probVec, 0, sizeof(REAL) * vocabLen); REAL colSum = 0;
			for (unsigned long col = 0; col < vocabLen; ++col)
			{
				REAL sum = 0;
				for (int i = 0; i < i_dim; ++i) sum += wordVec[row][i] * contextVec[col][i];
				probVec[col] = exp(sum + contextBias[col] + wordBias[row]); colSum += probVec[col];
			}
			for (unsigned long col = 0; col < vocabLen; ++col) trgtProbTransMat[row][col] = (probVec[col] / colSum);
		}

		ifsinfo1.seekg(0, std::ios::beg);
		for (int k = 0;k < num1 + 1;k++)
			getline(ifsinfo1, str_seqName);
		vec_10.push_back(str_seqName);
		vec_info.clear();
		split(str_seqName, ",", vec_info);
		alldiff = "";
		alldiff.append(vec_info[0]).append(",").append(vec_dis[vec_dis.size() - 1 - i]);


		for (int j = 0;j < curnum.size();j++)
		{

			stringstream ss;
			ss << curnum[curnum.size() - 1 - j];
			int num;
			ss >> num;
			ifsData.seekg(sizeof(REAL) * vocabLen * (2 * i_dim + 2) * num, std::ios::beg);

			memset(wordBias, 0, sizeof(REAL) * vocabLen); memset(contextBias, 0, sizeof(REAL) * vocabLen);
			for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx)
			{
				memset(wordVec[vocabIdx], 0, sizeof(REAL) * i_dim);
				memset(contextVec[vocabIdx], 0, sizeof(REAL) * i_dim);
			}
			for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx) ifsData.read((char*)wordVec[vocabIdx], sizeof(REAL)*i_dim);
			for (unsigned long vocabIdx = 0; vocabIdx < vocabLen; ++vocabIdx) ifsData.read((char*)contextVec[vocabIdx], sizeof(REAL)*i_dim);
			ifsData.read((char*)wordBias, sizeof(REAL)*vocabLen);
			ifsData.read((char*)contextBias, sizeof(REAL)*vocabLen);

			for (unsigned long row = 0; row < vocabLen; ++row)
			{
				memset(probVec, 0, sizeof(REAL) * vocabLen); REAL colSum = 0;
				for (unsigned long col = 0; col < vocabLen; ++col)
				{
					REAL sum = 0;
					for (int i = 0; i < i_dim; ++i) sum += wordVec[row][i] * contextVec[col][i];
					probVec[col] = exp(sum + contextBias[col] + wordBias[row]); colSum += probVec[col];
				}
				for (unsigned long col = 0; col < vocabLen; ++col) srcProbTransMat[row][col] = (probVec[col] / colSum);
			}





			REAL diff = 0.0;
			for (unsigned long row = 0; row < vocabLen; ++row)
			{
				for (unsigned long col = 0; col < vocabLen; ++col)
					diff += fabs(trgtProbTransMat[row][col] - srcProbTransMat[row][col]);
			}


			stringstream ss1;
			ss1 << diff;
			string currdiff;
			ss1 >> currdiff;
			alldiff.append(",").append(currdiff);
			Mat_dist->set(i + 1, j + 1, diff);
		}
		outMatFile << alldiff << endl;
	}



	for (int r = 0;r < 11;r++) {
		for (int c = 0;c < 11;c++) {
			stringstream ss2;
			ss2 << Mat_dist->get(r, c);
			jjj = "";
			ss2 >> jjj;
			cout << Mat_dist->get(r, c) << ",";
		}
		cout << endl;
	}

	Mat_pcoa = Mat_dist->MDS_UCF(2, 60);


	for (int r = 0;r < 11;r++) {
		for (int c = 0;c < 2;c++) {
			cout << Mat_pcoa->get(r, c) - Mat_pcoa->get(0, c) << ",";
		}
		cout << endl;
	}



	for (int r = 0;r < 10;r++) {
		outinfoall = vec_dis[9 - r];
		for (int c = 0;c < 2;c++) {
			stringstream ss2;
			ss2 << Mat_pcoa->get(10 - r, c) - Mat_pcoa->get(0, c);
			jjj = "";
			ss2 >> jjj;
			outinfoall.append(",").append(jjj);

		}
		outinfoall.append(",").append(vec_10[r]);
		outinfoFile << outinfoall << endl;
	}
	outMatFile.close();
	outinfoFile.close();
	curnum.clear();

	endTime = clock();
	std::cout << "Query Time Elapsed: " << ((float)endTime - (float)startTime) / CLOCKS_PER_SEC << " seconds" << std::endl;//
	for (unsigned long lineCnt = 0; lineCnt < vocabLen; ++lineCnt) { delete[] wordVec[lineCnt]; delete[] contextVec[lineCnt]; delete[] srcProbTransMat[lineCnt]; delete[] trgtProbTransMat[lineCnt]; }
	delete[] probVec; delete[] wordVec; delete[] contextVec; delete[] wordBias; delete[] contextBias; delete[] srcProbTransMat; delete[] trgtProbTransMat;


	return 0;
}



