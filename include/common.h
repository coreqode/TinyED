#ifndef COMMON_H
#define COMMON_H

#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace TinyED
{
    typedef double Real;

    template<typename T>
    using MatSparse = Eigen::SparseMatrix<T>;

    template<unsigned int N, unsigned int M>
    using MatS = Eigen::Matrix<Real, N, M, Eigen::DontAlign>;

    template<unsigned int N>
    using VecS = Eigen::Matrix<Real, N, 1, Eigen::DontAlign>;

    using MatD = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

    using VecD = Eigen::Matrix<Real, Eigen::Dynamic, 1>;

}

#endif //COMMON_H
