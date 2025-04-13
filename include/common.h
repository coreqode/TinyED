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

    template<unsigned int N, unsigned int M>
    using Mat = Eigen::Matrix<Real, N, M, Eigen::DontAlign>;

    template<unsigned int N>
    using VecS = Eigen::Matrix<Real, N, 1, Eigen::DontAlign>;

    template<unsigned int N>
    using Vec = Eigen::Matrix<Real, N, 1, Eigen::DontAlign>;


    using Mat2 = Eigen::Matrix<Real, 2, 2, Eigen::DontAlign>;
    using Mat3 = Eigen::Matrix<Real, 3, 3, Eigen::DontAlign>;
    using Mat32 = Eigen::Matrix<Real, 3, 2, Eigen::DontAlign>;
    using Mat23 = Eigen::Matrix<Real, 2, 3, Eigen::DontAlign>;
    using Vec3 = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
    using Vec2 = Eigen::Matrix<Real, 2, 1, Eigen::DontAlign>;
    using VecN = Eigen::Matrix<Real, Eigen::Dynamic, 1, Eigen::DontAlign>;

    using MatN = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::DontAlign>;
    using Mat3N = Eigen::Matrix<Real, 3, Eigen::Dynamic, Eigen::DontAlign>;
    using Mat2N = Eigen::Matrix<Real, 2, Eigen::Dynamic, Eigen::DontAlign>;


    template <class T>
    using Map = Eigen::Map<T>;

    template <class T>
    inline Map<const VecN> FlattenMatrixN3(const T& vector)
    {
        return Map<const VecN>(vector.data(), 3 * vector.cols());
    }


    template <class T>
    inline Map<VecN> FlattenMatrixN3(T& vector)
    {
        return Map<VecN>(vector.data(), 3 * vector.cols());
    }


    template <class T>
    inline Map<const Mat3N> GroupMatrixN3(const T& vector)
    {
        return Map<const Mat3N>(vector.data(), 3, vector.size() / 3);
    }


}

#endif //COMMON_H
