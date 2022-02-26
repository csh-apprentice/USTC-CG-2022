#pragma once
// todo 5: change the class in to a template class
#ifndef _MATRIX_SPARSE_H
#define _MATRIX_SPARSE_H
#include<iostream>
#include <cstdio>
#include <vector>
#include <assert.h>
#include <fstream>
//#include "matrix_sparse.h"

#define my_epsilon 1e-5
using namespace std;
template<class T>
class matrix_sparse;

template<class T>
class matrix_base {
public:
    friend class matrix_sparse<T>;
    // default constructor
    // https://en.cppreference.com/w/cpp/language/constructor
    matrix_base()// 
    {/*.todo 1..*/
     //default: Initialize a 4*4 Identity Matrix 
        rows = cols = 4;
        data = new T[rows * cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                if (j == i) data[i * cols + j] = T(1);
                else data[i * cols + j] = T(0);
            }
    }

    // constructor with initilizer list
    matrix_base(int rows, int cols) 
    {/*.todo 1..*/
     //default: Initialize a Zero matrix with rows and cols provided.
        this->rows = rows;
        this->cols = cols;
        data = new T[rows * cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[i*cols+j] = T(0);
    }
    //read from file
    matrix_base(const std::string& file)
    {
        assert(Read_from_file(file));
    }

    //read from list, if data_vec.size()>rows*cols, we intercept the data, otherwise we fill the matrix with zero
    matrix_base(int rows,int cols,vector<T> data_vec)
    {
        assert(rows > 0 && cols > 0);
        int index = 0;    //the index of the entry of the matrix
        data = new T[rows * cols];
        this->rows = rows;
        this->cols = cols;
        //cout << "safe here" << endl;
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {
                data[index] = data_vec.size() ? data_vec[index] : T(0);
                index++;
            }             
    }
    //copy initializer
    matrix_base(const matrix_base& y) 
    {
        rows = y.nrow();
        cols = y.ncol();
        data = new T[rows * cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                data[i * cols + j] = y.data[i*cols+j];
    
    }
    //initilize from a sparse matrix
    matrix_base(matrix_sparse<T> s)
    {
        rows = s.nrow();
        cols = s.ncol();
        data = new T[rows * cols];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                if (s.s_data.find(make_pair(i, j)) == s.s_data.end())
                    data[i * cols + j] = T(0);
                else data[i * cols + j] = s.s_data[make_pair(i, j)];

    }

    // desctructor
    // https://en.cppreference.com/w/cpp/language/destructor
    ~matrix_base() 
    {/*.todo 2..*/
        delete[]data;
    }

    int nrow() const { return rows; }
    int ncol() const { return cols; }

    // operator overloding
   T& operator()(int r, int c)   
    { /* todo 3: particular entry of the matrix*/ 
      //try to return thr enttry[i,j] of the matrix, we shall first judge if the index is out of range, here we use assert;
       assert(!(r < 0 || r >= this->rows || c < 0 || c >= this->cols));
        return this->data[r * this->cols + c];
    }
    T& operator[](int index) 
    { 
        /* todo 3: particular entry of the matrix*/
        assert(!(index < 0 && index >= this->rows * this->cols));
            
        return this->data[index];
    }


    matrix_base<T> col(int c) {
       /* todo 4: particular column of the matrix*/ 
        assert(c >= 0 && c < this->ncol());
        vector<T> data_vec;
        for (int i = 0; i < this->rows; i++)
        {
            data_vec.push_back((*this)(i, c));
            cout << (*this)(i, c) << endl;

        }
        matrix_base colmatrix = matrix_base(this->rows, 1, data_vec);
        //colmatrix.print();
        return colmatrix;
    
    }
    matrix_base row(int r) 
    {
        /* todo 4: particular row of the matrix*/
        assert(r >= 0 && r < this->nrow());
        vector<T> data_vec;
        for (int i = 0; i < this->cols; i++)
        {
            data_vec.push_back((*this)(r, i));        
        }
        matrix_base rowmatrix = matrix_base(1, this->cols, data_vec);
        rowmatrix.print();
        return rowmatrix; 
    }
    matrix_base submat(int rl, int rh, int cl, int ch) const 
    { 
        /* todo 4: return a sub-matrix specified by the input parameters*/
        assert(rh >= rl && ch >= cl);
        assert(rl >= 0 && rh < this->nrow());
        assert(cl >= 0 && ch < this->ncol());
        vector<T> data_sub;
        for (int i = rl; i <= rh; i++)
            for (int j = cl; j <= ch; j++)
                data_sub.push_back(this->data[i * this->cols + j]);
        int sub_row = rh - rl + 1, sub_col = ch - cl + 1;
        return matrix_base(sub_row,sub_col,data_sub); 
    }

    // constant alias
    matrix_base& operator= (const matrix_base& rhs) 
    {/*.todo 3..*/
        this->cols = rhs.ncol();
        this->rows = rhs.nrow();
        //this->data = new T[this->nrow() * this->ncol()];
        for (int i = 0; i < this->nrow(); i++)
            for (int j = 0; j < this->ncol(); j++)
                this->data[i * this->ncol() + j] = rhs.data[i * this->ncol() + j];
        return *this;
    }
    
    
    matrix_base operator+ (const matrix_base& rhs)
    {/*.todo 3..*/
        assert(this->rows == rhs.nrow() && this->cols == rhs.ncol());
        vector<T> sum_data;
        for (int i = 0; i < this->nrow(); i++)
            for (int j = 0; j < this->ncol(); j++)
            {
                sum_data.push_back(this->data[i * this->ncol() + j] + rhs.data[i * this->ncol() + j]);
            }
        return matrix_base(rows, cols, sum_data);
    }
    
    matrix_base operator- (const matrix_base& rhs) 
    {/*.todo 3..*/
        assert(this->rows == rhs.nrow() && this->cols == rhs.ncol());
        vector<T> sub_data;
        for (int i = 0; i < this->nrow(); i++)
            for (int j = 0; j < this->ncol(); j++)
            {
                sub_data.push_back(this->data[i * this->ncol() + j] - rhs.data[i * this->ncol() + j]);
            }
        return matrix_base(rows, cols, sub_data);
     }
    matrix_base operator* (const matrix_base& rhs) 
    {/*.todo 3..*/
        assert(this->ncol() == rhs.nrow());
        vector<T> mul_data;
        mul_data.resize(rhs.ncol() * this->nrow());
        for (int i = 0; i < this->nrow(); i++)
            for(int j=0;j< rhs.ncol();j++)
            {
                T num = T(0);
                for (int k = 0; k < this->ncol(); k++)        
                    num += this->data[i * cols + k] * rhs.data[k * rhs.ncol() + j];
                //cout << "index is" << i * cols + j << endl;
                mul_data[i * rhs.ncol() + j] = num;
                
            }
        return matrix_base(this->nrow(),rhs.ncol(),mul_data);
    }
    matrix_base inverse()
    {
        assert(rows==cols);
        int size = rows;
        int i, j, k;
        matrix_base m(*this);
        vector<T> s;
        s.resize(m.nrow() * m.ncol());
        for (int index = 0; index < m.nrow() * m.ncol();index+=size+1)
            s[index] = T(1);
        for (i = 0; i < size-1; i++) {//find the max entry in thr lower triangle
            int pivot = i;
            T pivotsize = m.data[i*size+i];
            if (pivotsize < 0)pivotsize = -pivotsize;
            for (j = i + 1; j < size; j++) {
                T tmp = m.data[j*size+i];
                if (tmp < 0)tmp = -tmp;
                if (tmp > pivotsize) {
                    pivot = j;
                    pivotsize = tmp;
                }
            }
            if (pivotsize == 0) {
                throw - 1;
            }
            if (pivot != i)//swap two rows, enable the entry positioned in diagnal position to be the maximize entry among the column
            {
                for (j = 0; j < size; j++) {
                    T tmp;
                    tmp = m.data[i*size+j];
                    m.data[i*size+j] = m.data[pivot*size+j];
                    m.data[pivot*size+j] = tmp;
                    tmp = s[i*size+j];
                    s[i*size+j] = s[pivot*size+j];
                    s[pivot*size+j] = tmp;
                }
            }
            //set the entries located in the lower triangle = 0
            for (j = i + 1; j < size; j++) {
                T f = m.data[j*size+i] / m.data[i*size+i];
                for (k = 0; k < size; k++) {
                    m.data[j*size+k] -= m.data[i*size+k] * f;
                    s[j*size+k] -= s[i*size+k] * f;
                }
            }
        }

        for (i = size-1; i >= 0; i--) {
            T f;
            f = m.data[i*size+i];
            if (fabs(f) <my_epsilon) { throw - 1; }
            for (j = 0; j < size; j++) {//set the entry located in diagnal position = 1
                m.data[i*size+j] /= f;
                s[i*size+j] /= f;
            }
            for (j = 0; j < i; j++) {
                f = m.data[j*size+i];
                for (k = 0; k < size; k++) {
                    m.data[j*size+k] -= f * m.data[i*size+k];
                    s[j*size+k] -= f * s[i*size+k];
                }
            }
        }       
        return matrix_base(size,size,s);
    }
    
    
    
    matrix_base operator/ (const matrix_base& rhs) {
        /*.todo 3..*/ 
        matrix_base cop = rhs; //first copy the matrix
        return (*this)*cop.inverse();
    }

    matrix_base operator+= (const matrix_base& rhs) 
    {/*.todo 3..*/ 
        (*this) = (*this) + rhs;
        return (*this); 
    }
    matrix_base operator-= (const matrix_base& rhs) 
    {/*.todo 3..*/ 
        (*this) = (*this) - rhs;
        return (*this); 
    }
    matrix_base operator*= (const matrix_base& rhs) 
    {/*.todo 3..*/ 
        (*this) = (*this) * rhs;
        return (*this);
    }
    matrix_base operator/= (const matrix_base& rhs) 
    {/*.todo 3..*/ 
        (*this) = (*this) / rhs;
        return (*this);
    }

    matrix_base operator+ (T v) 
    {/*.todo 3..*/
        matrix_base cop = (*this);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                cop.data[i * cols + j] += v;
        return cop;
    }
    matrix_base operator- (T v) 
    {/*.todo 3..*/
        matrix_base cop = (*this);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                cop.data[i * cols + j] -= v;
        return cop;
    }
    matrix_base operator* (T v) 
    {/*.todo 3..*/
        matrix_base cop = (*this);
        for (int i = 0; i < rows; i++)
            for (int j = cols;j<cols; j++)
                cop.data[i * cols + j] *= v;
        return cop; 
    }
    matrix_base operator/ (T v) 
    {/*.todo 3..*/
        assert(fabs(v) > my_epsilon);
        matrix_base cop = (*this);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                cop.data[i * cols + j] /= v;
        return cop;
    }

    matrix_base operator+= (T v) 
    {/*.todo 3..*/
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                this->data[i * cols + j] += v;
        return (*this);
    }
    matrix_base operator-= (T v) 
    {/*.todo 3..*/
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                this->data[i * cols + j] -= v;
        return (*this);
    }
    matrix_base operator*= (T v) 
    {/*.todo 3..*/
        for (int i = 0; i < rows; i++)
            for (int j = 0;j<cols; j++)
                this->data[i * cols + j] *= v;
        return (*this);
    }
    matrix_base operator/= (T v) 
    {/*.todo 3..*/
        assert(fabs(v) > my_epsilon);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                this->data[i * cols + j] /= v;
        return (*this);
    }

    friend ostream& operator<<(ostream& out, T& t)    //输出流重载声明及实现
    {
        return out << "data   is   " << t;
    } //--------------------------------------------
    friend istream& operator>>(istream& in, T& t)      //输入流重载声明及实现
    {
        return in >> t;
    }//---------------------------------------------

    void print() const {
        printf("this matrix has size (%d x %d)\n", rows, cols);
        printf("the entries are:\n");
        /* todo 4: print all the entries of the matrix */
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
                cout << data[i * cols + j] << " ";
            cout << endl;
        }
    }
    bool Read_from_file(const std::string& file)
    {
    
        std::ifstream inp;
        inp.open(file.c_str());
        if (!inp.is_open()) {
            std::cout << "ERROR::matrix_base::ReadFromFile:" << std::endl
                << "\t" << "file [" << file << "] opens failed" << std::endl;
            return false;
        }
        inp >> rows >> cols;            //read the dimension of the matrix
        assert(rows >= 0 && cols >= 0);
        data = new T[rows * cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
            {   
                inp>>data[i * cols + j];           //read the element to constuct the initial matrix
            }

        inp.close();
        return true;
    
    }
    
    
private:
    int rows=-1, cols=-1;
    T* data;
};
#endif