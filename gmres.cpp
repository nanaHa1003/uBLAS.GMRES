#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

template<class type>
inline static void rotate(type c, type s, type &r, type &h)
{
    type tmp = c * r + s * h;
    h        = c * h - s * r;
    r        = tmp;
}

template<class type>
int gmres(ublas::vector<type> &y, ublas::matrix<type> &A, ublas::vector<type> &x, type tol, size_t m)
{
    using namespace ublas;

    typedef vector<type> Vec;
    typedef matrix<type> Mat;

    int size = y.size();
    unit_vector<type> e_1(size + 1, 0);

    m = (m > size) ?size :m;

    // Lower case name for vector
    Vec x_i(x);
    Vec r_i = y - prod(A, x);
    Vec v_i = r_i / norm_2(r_i);
    Vec e_i = e_1 * norm_2(r_i);

    // Givens rotation args
    Vec c(m + 1, 0);
    Vec s(m + 1, 0);

    // Upper case name for matrix, V = V^t, R = Q^t * H
    Mat V(m, size, 0);
    Mat R(size + 1, m, 0);

    for(int i = 0; i < size; i++)
    {
        type beta = norm_2(r_i);
        Vec  v_i  = r_i / beta;
        Vec  e_i  = e_1 * beta;

        int dim = 0;
        V.clear();
        for(int j = 0; j < m; j++)
        {
            row(V, j) = v_i;
            v_i       = prod(A, v_i);
            for(int k = 0; k <= j; k++)
            {
                R(k, j) = inner_prod(v_i, row(V, k));
                v_i    -= R(k, j) * row(V, k);
            }

            // Re-orthogonalization
            #pragma omp parallel for
            for(int k = 0; k <= j; k++)
            {
                type tmp   = inner_prod(v_i, row(V, k));
                row(V, k) -= tmp * v_i;
            }

            R(j + 1, j) = norm_2(v_i);

            dim++;
            if(R(j + 1, j) > tol)
            {
                v_i /= R(j + 1, j);
            }
            else
            {
                R(j + 1, j) = (type) 0;
                v_i.clear();
                break;
            }
        }

        // Apply givens rotation
        for(int j = 0; j < m; j++)
        {
            type r = std::sqrt(R(j, j) * R(j, j) + R(j + 1, j) * R(j + 1, j));
            if(r > tol)
            {
                c(j)   = R(j, j)     / r;
                s(j)   = R(j + 1, j) / r;
            }
            else
            {
                c(j) = (type) 1;
                s(j) = (type) 0;
            }

            rotate(c(j), s(j), R(j, j), R(j + 1, j));
            rotate(c(j), s(j), e_i(j), e_i(j + 1));
        }

        // Solve for y_i
        Vec y_i(solve(subrange(R, 0, dim, 0, dim), subrange(e_i, 0, dim), upper_tag()));

        // Update x
        x_i += prod(y_i, subrange(V, 0, dim, 0, size));
        r_i  = y - prod(A, x_i);

        if(norm_2(r_i) < tol) break;
    }

    x = x_i;
    return 0;
}

int main(int argc, char **argv)
{
    std::cout.precision(15);

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
    if(gmres(y, A, x, 1e-12, y.size()) != 0)
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

