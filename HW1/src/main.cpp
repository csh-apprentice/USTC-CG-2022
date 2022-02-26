#include "matrix_base.h"
#include "matrix_sparse.h"
#include <iostream>
#include <chrono>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;
using namespace chrono;

//typedef Eigen::Triplet<double> T;

int main() {
    int size = 4;
    cout << "----------------------Followings are the display of the dense matrix:--------------------" << endl;

    matrix_base<double> A(size,size), B(size,size-1);


    // todo 6: fill A anb B with random numbers
    for(int i=0; i<A.nrow(); i++)
        for (int j = 0; j < A.ncol(); j++)
            A(i, j) = rand()/10;
       
    for (int i = 0; i < B.nrow(); i++)
        for (int j = 0; j < B.ncol(); j++)
            B(i, j) = rand()/10;

    // todo 7: benchmark with runtime, using std::chrono
    // https://en.cppreference.com/w/cpp/chrono
    steady_clock::time_point t1 = steady_clock::now();

    matrix_base<double> C = A*B;
    cout << "Here is the Matrix A in my own implemention: " << endl;
    A.print();
    cout << "Here is the Matrix B in my own implemention: " << endl;
    B.print();
    cout << "Here is the Matrix C in my own implemention: " << endl;
    C.print();
    steady_clock::time_point t2 = steady_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "It took me " << time_span.count() << " seconds in my own implemented method." << endl;
    // todo 8: use Eigen and compare

    MatrixXd eigen_A = Eigen::MatrixXd::Identity(size, size);
    MatrixXd eigen_B = Eigen::MatrixXd::Zero(size, size-1);
    for (int i = 0; i < A.nrow(); i++)
        for (int j = 0; j < A.ncol(); j++)
            eigen_A(i, j) = A(i,j);

    for (int i = 0; i < B.nrow(); i++)
        for (int j = 0; j < B.ncol(); j++)
            eigen_B(i, j) = B(i, j);

    steady_clock::time_point t3 = steady_clock::now();
    MatrixXd eigen_C = eigen_A * eigen_B;
    cout << "Here is the Matrix A in Eigen: " << endl;
    cout << eigen_A << endl;
    cout << "Here is the Matrix B in Eigen: " << endl;
    cout << eigen_B << endl;
    cout << "Here is the Matrix C in Eigen: " << endl;
    cout << eigen_C << endl;
    steady_clock::time_point t4 = steady_clock::now();
    time_span = duration_cast<duration<double>>(t4 - t3);
    std::cout << "It took me " << time_span.count() << " seconds in Eigen Library." << endl;

    cout << "----------------------Followings are the display of the sparse matrix:--------------------" << endl;
    matrix_sparse<double> sparse_A(A),sparse_B(size),sparse_C(size,size);
    t1 = steady_clock::now();

    try {
        sparse_C = sparse_B / sparse_A;
    }
    catch (int)
    {
        cerr << "the matrix doesn't have an inverse version!" << endl;
    }
    
    sparse_A.print();
    sparse_B.print();
    sparse_C.print();
    t2 = steady_clock::now();
    std::cout << "It took me " << time_span.count() << " seconds in my own implemented method." << endl;
  
   

    return 0;
}