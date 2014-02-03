#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

template<class type>
int gmres(ublas::vector<type> &y, ublas::matrix<type> &A, ublas::vector<type> &x)
{
    using namespace ublas;

    typedef vector<type>                Vec;
    typedef matrix<type>                Mat;
    typedef matrix_range <matrix<type>> SubMat;

    int size = y.size();
    unit_vector<type> e_1(size + 1, 0);

    // Lower case name for vector
    Vec x_i(x);
    Vec r_i = y - prod(A, x);
    Vec v_i = r_i / norm_2(r_i);
    Vec e_i = e_1;

    // Givens rotation args
    Vec c(size + 1, 0);
    Vec s(size + 1, 0);

    // Upper case name for matrix, V = V^t, R = Q^t * H
    Mat V(size, size, 0);
    Mat R(size + 1, size, 0);

    // size can be subsitude by m
    for(int i = 0; i < size; i++)
    {
        type beta = norm_2(r_i);
        std::cout << "beta = " << beta << std::endl;
        if(e_i(i) / norm_2(y) < 1e-8) break;

        // Add v_i to V
        row(V, i) = v_i;

        // Build R_i
        SubMat V_i(V, range(0, i + 1), range(0, size));
        SubMat R_i(R, range(0, i + 2), range(0, i + 1));
        Vec    tmp1(prod(A, v_i));
        Vec    tmp2(prod(V_i, tmp1));

        #pragma omp parallel for
        for(int j = 0; j < tmp2.size(); j++)
        {
            R(j, i) = tmp2(j);
        }

        // Update v_i
        noalias(v_i) = prod(A, row(V, i));
        for(int j = 0; j <= i; j++)
        {
            v_i -= R(j, i) * row(V, i);
        }

        R(i + 1, i) = norm_2(v_i);
        v_i        /= R(i + 1, i);

        // Apply Givens rotation
        for(int j = 0; j < i; j++)
        {
            type temp   = c(j) * R(j, i) - s(j) * R(j + 1, i);
            R(j + 1, i) = c(j) * R(j + 1, i) + s(j) * R(j ,i);
            R(j, i)     = temp;
        }

        type ratio = std::sqrt(R(i + 1, i) * R(i + 1, i) + R(i, i) * R(i, i));

        c(i) =  R(i, i) / ratio;
        s(i) = -R(i + 1, i) / ratio;

        // R(i, i) = sqrt(r^2 + h^2)
        R(i, i)     = ratio;
        R(i + 1, i) = (type) 0;

        // Apply rotation on e_i
        type rot_tmp = c(i) * e_i(i) - s(i) * e_i(i + 1);
        e_i(i + 1)   = c(i) * e_i(i + 1) + s(i) * e_i(i);
        e_i(i)       = rot_tmp;

        // Solve for y_i
        Vec y_i(solve(subrange(R, 0, i + 1, 0, i + 1), beta * subrange(e_i, 0, i + 1), upper_tag()));

        // Update x
        noalias(x_i) += prod(y_i, V_i);
        noalias(r_i)  = y - prod(A, x_i);
    }

    x = x_i;
    return 0;
}

int main(int argc, char **argv)
{
    int size, nnz;
    std::string filename;
    if(argc == 2) filename = argv[1];
    else filename = "prob.mtx";

    std::ifstream prob(filename);
    if(!prob.is_open()) exit(0);
    prob >> size;
    prob >> size >> nnz;

    std::cout << "Problem size : " << size << std::endl;

    ublas::vector<double> x(size, 0.0);
    ublas::vector<double> y(size, 1.0);
    ublas::matrix<double> A(size, size, 0.0);

    // Set up problem
    for(int i = 0; i < nnz; i++)
    {
        int tmpRow, tmpCol;
        double tmpVal;

        prob >> tmpRow >> tmpCol >> tmpVal;
        A(tmpRow - 1, tmpCol - 1) = tmpVal;
    }

    // GMRES
    if(gmres(y, A, x) != 0)
    {
        std::cout << "GMRES not converged" << std::endl;
    }
    else
    {
        double r = norm_2(y - prod(A, x));
        std::cout << "r   = " << r << std::endl;
        std::cout << "err = " << r / norm_2(y) << std::endl;
    }

    std::cout << x << std::endl;

    return 0;
}

