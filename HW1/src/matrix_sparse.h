#pragma once
#ifndef _MATRIX_BASE_H
#define _MATRIX_BASE_H
#include "matrix_base.h"
#include <iostream>
#include <map>
#include <vector>
#include <utility>


template<class T>
class matrix_base;

template<class T>
class matrix_sparse
{
public:
	friend class matrix_base<T>;
	
	//default initilaizer, create a Zero 4*4 sparse matrix
	matrix_sparse() 
	{
		s_rows = s_cols = 4;
	}

	//create a Identity matrix
	matrix_sparse(int size)
	{
		assert(size > 0);
		s_rows = s_cols = nozero= size;
		for(int i=0;i<size;i++)
			s_data.insert(map<pair<int, int>, T>::value_type(make_pair(i, i), T(1)));
	
	}

	//initilize a Zero matrix with given row and col
	matrix_sparse(int r, int c)
	{
		s_rows = r;
		s_cols = c;
	}

	//initialize a sparse marix with all the entry is v
	matrix_sparse(int r, int c, T v)
	{
		s_rows = r;
		s_cols = c;
		nozero = r * c;
		for (int i = 0; i < r; i++)
			for (int j = 0; j < c; j++)
				s_data.insert(map<pair<int, int>, T>::value_type(make_pair(i, j), v));
	}

	matrix_sparse(int r, int c, map<pair<int, int>, T> data_2)
	{
		this->s_rows = r;
		this->s_cols = c;
		this-> nozero = data_2.size();
		
		copy(data_2.begin(), data_2.end(), inserter(s_data, s_data.begin()));
		
	
	}

	//initiaize a Zero matrix with a dense matrix from the matrix_base
	matrix_sparse(matrix_base<T> dense)
	{
		s_rows = dense.nrow();
		s_cols = dense.ncol();
		for(int i=0;i<s_rows;i++)
			for (int j = 0; j < s_cols; j++)
			{
				if (fabs(dense(i, j)) > my_epsilon)
				{
					s_data.insert(map<pair<int, int>, T>::value_type(make_pair(i, j), dense(i, j)));
					nozero++;
				}
			}
	}
	

	~matrix_sparse() {}
	int nrow() const { return s_rows; }
	int ncol() const { return s_cols; }
	int get_nozero() const { return nozero; }

	// operator overloding
	T operator()(int r, int c)
	{ /* todo 3: particular entry of the matrix*/
	  //try to return thr enttry[i,j] of the matrix, we shall first judge if the index is out of range, here we use assert;
		assert(!(r < 0 || r >= this->s_rows || c < 0 || c >= this->s_cols));
		if (s_data.find(make_pair(r, c)) != s_data.end())
			return this->s_data[make_pair(r, c)];
		else return T(0);
	}
	T operator[](int index)
	{
		/* todo 3: particular entry of the matrix*/
		assert(!(index < 0 && index >= this->s_rows * this->s_cols));
		int r = index / s_cols;
		int c = index - r * s_cols;
		return (*this)(r,c);
	}

	matrix_sparse<T> col(int c)
	{
		assert(c >= 0 && c < this->ncol());
		map<pair<int, int>, T> data_2;
		for (auto iter = s_data.begin(); iter != s_data.end(); iter++)
		{
			//cout << (iter->first).second <<"ssss"<< endl;
			if ((iter->first).second == c)
			{
				data_2.insert(map<pair<int, int>, T>::value_type(make_pair((iter->first).first, 0), iter->second));
			}
		}
		
		return matrix_sparse(s_rows, 1, data_2);
	}

	matrix_sparse<T> row(int r)
	{
		assert(r >= 0 && r < this->nrow());
		map<pair<int, int>, T> data_2;
		for(auto iter=s_data.begin();iter!=s_data.end();iter++)
			if((iter->first).first==r)
				data_2.insert(map<pair<int, int>, T>::value_type(make_pair(0, (iter->first).second), iter->second));
		return matrix_sparse(1, s_cols, data_2);
	}

	matrix_sparse submat(int rl, int rh, int cl, int ch)
	{
		/* todo 4: return a sub-matrix specified by the input parameters*/
		assert(rh >= rl && ch >= cl);
		assert(rl >= 0 && rh < this->nrow());
		assert(cl >= 0 && ch < this->ncol());
		map<pair<int, int>, T> data_2;
		for (auto iter = s_data.begin(); iter != s_data.end(); iter++)
			if ((iter->first).first >= rl && (iter->first).first <= rh && (iter->first).second >= cl && (iter->first).second <= ch)
			{
				//data_2.insert(map<pair<int, int>, T>::value_type(make_pair((iter->first).first- rl , (iter->first).second) - cl), iter->second));
				data_2.insert(map<pair<int, int>, T>::value_type(make_pair((iter->first).first - rl, (iter->first).second-cl), iter->second));
			}
		int sub_row = rh - rl + 1, sub_col = ch - cl + 1;
	    return matrix_sparse(sub_row, sub_col, data_2);
		
	}
	// constant alias
	matrix_sparse& operator= (const matrix_sparse& rhs)
	{/*.todo 3..*/
		this->s_cols = rhs.ncol();
		this->s_rows = rhs.nrow();
		this->nozero = rhs.nozero;
		this->s_data = rhs.s_data;
		
		return *this;
	}
	matrix_sparse operator+ (const matrix_sparse& rhs)
	{/*.todo 3..*/
	    assert(this->s_rows == rhs.nrow() && this->s_cols == rhs.ncol());
		map<pair<int, int>, T> sum_data;
		auto cop = rhs.s_data;
		for(int i=0;i<s_rows;i++)
			for (int j = 0; j < s_cols; j++)
			{
				T num1 = 0, num2 = 0;
				if (s_data.find(make_pair(i, j)) != s_data.end()) num1 = s_data[make_pair(i, j)];
				if (rhs.s_data.find(make_pair(i, j)) != rhs.s_data.end())  num2 = cop[make_pair(i, j)];
				if (fabs(num1 + num2) > my_epsilon) sum_data.insert(map<pair<int, int>, T>::value_type(make_pair(i, j), num1 + num2));
				
			}
		return matrix_sparse(s_rows, s_cols, sum_data);
	}

	matrix_sparse operator- (const matrix_sparse& rhs)
	{/*.todo 3..*/
		assert(this->s_rows == rhs.nrow() && this->s_cols == rhs.ncol());
		map<pair<int, int>, T> sum_data;
		auto cop = rhs.s_data;
		for (int i = 0; i < s_rows; i++)
			for (int j = 0; j < s_cols; j++)
			{
				T num1 = 0, num2 = 0;
				if (s_data.find(make_pair(i, j)) != s_data.end()) num1 = s_data[make_pair(i, j)];
				if (rhs.s_data.find(make_pair(i, j)) != rhs.s_data.end())  num2 = cop[make_pair(i, j)];
				if (fabs(num1 - num2) > my_epsilon) sum_data.insert(map<pair<int, int>, T>::value_type(make_pair(i, j), num1 - num2));

			}
		return matrix_sparse(s_rows, s_cols, sum_data);
	}

	matrix_sparse operator* (const matrix_sparse& rhs)
	{/*.todo 3..*/
		assert(this->ncol() == rhs.nrow());
		matrix_base<T> dense1 = (*this),dense2=rhs;
		return matrix_sparse(dense1*dense2);
	}

	matrix_sparse operator/ (const matrix_sparse& rhs) {
		/*.todo 3..*/
		matrix_base<T>  dense1=(*this),dense2 = rhs; //first copy the matrix
		return matrix_sparse(dense1 * dense2.inverse());
	}

	matrix_sparse operator+= (const matrix_sparse& rhs)
	{/*.todo 3..*/
		(*this) = (*this) + rhs;
		return (*this);
	}
	matrix_sparse operator-= (const matrix_sparse& rhs)
	{/*.todo 3..*/
		(*this) = (*this) - rhs;
		return (*this);
	}
	matrix_sparse operator*= (const matrix_sparse& rhs)
	{/*.todo 3..*/
		(*this) = (*this) * rhs;
		return (*this);
	}

	matrix_sparse operator/= (const matrix_sparse& rhs)
	{/*.todo 3..*/
		(*this) = (*this) / rhs;
		return (*this);
	}

	matrix_sparse operator+ (T v)
	{/*.todo 3..*/
		matrix_sparse cop = (*this);
		return cop+matrix_sparse(s_rows,s_cols,v);
	}

	matrix_sparse operator- (T v)
	{/*.todo 3..*/
		matrix_sparse cop = (*this);
		return cop - matrix_sparse(s_rows, s_cols, v);
	}

	matrix_sparse operator* (T v)
	{/*.todo 3..*/
		auto data2 = s_data;
		for (auto iter = data2.begin(); iter != data2.end(); iter++)
			data2[iter->first] *= v;
		return matrix_sparse(s_rows, s_cols, data2);
	}

	matrix_sparse operator/ (T v)
	{/*.todo 3..*/
		auto data2 = s_data;
		for (auto iter = data2.begin(); iter != data2.end(); iter++)
			data2[iter->first] /= v;
		return matrix_sparse(s_rows, s_cols, data2);
	}

	matrix_sparse operator+= (T v)
	{/*.todo 3..*/
		(*this) = (*this) + v;
		return (*this);
	}

	matrix_sparse operator-= (T v)
	{/*.todo 3..*/
		(*this) = (*this) - v;
		return (*this);
	}

	matrix_sparse operator*= (T v)
	{/*.todo 3..*/
		(*this) = (*this) * v;
		return (*this);
	}

	matrix_sparse operator/= (T v)
	{/*.todo 3..*/
		(*this) = (*this) / v;
		return (*this);
	}

	void print() const {
		printf("this matrix has size (%d x %d)\n", s_rows, s_cols);
		printf("the entries are:\n");
		/* todo 4: print all the entries of the matrix */
		for (int i = 0; i < s_rows; i++)
		{
			for (int j = 0; j < s_cols; j++)
			{
				auto iter = s_data.find(make_pair(i, j));
				if (iter == s_data.end())
					cout << T(0) << " ";
				else
					cout << iter->second << " ";
			} 
				
			cout << endl;
		}
	}
private:
	int s_rows = -1, s_cols = -1,nozero=0;
	map<pair<int, int>, T> s_data;   //we use map stucture to storge and do matrix calculation
};


#endif