#ifndef COMMON_H
#define COMMON_H

#include <iostream>
#include "Eigen/Dense"
#include "Eigen/Sparse"
typedef double Real;

using Eigen::VectorXd;
using Eigen::MatrixXd;

using Vector2r = Eigen::Matrix<Real, 2, 1, Eigen::DontAlign>;
using Vector3r = Eigen::Matrix<Real, 3, 1, Eigen::DontAlign>;
using Vector4r = Eigen::Matrix<Real, 4, 1, Eigen::DontAlign>;
using Vector5r = Eigen::Matrix<Real, 5, 1, Eigen::DontAlign>;
using Vector6r = Eigen::Matrix<Real, 6, 1, Eigen::DontAlign>;
using Vector6i = Eigen::Matrix<int, 6, 1, Eigen::DontAlign>;
using Vector18r = Eigen::Matrix<Real, 18, 1, Eigen::DontAlign>;
using Matrix2r = Eigen::Matrix<Real, 2, 2, Eigen::DontAlign>;
using Matrix3r = Eigen::Matrix<Real, 3, 3, Eigen::DontAlign>;
using Matrix4r = Eigen::Matrix<Real, 4, 4, Eigen::DontAlign>;
using Vector2i = Eigen::Matrix<int, 2, 1, Eigen::DontAlign>;
using AlignedBox2r = Eigen::AlignedBox<Real, 2>;
using AlignedBox3r = Eigen::AlignedBox<Real, 3>;
using AngleAxisr = Eigen::AngleAxis<Real>;
using Quaternionr = Eigen::Quaternion<Real, Eigen::DontAlign>;

using RowVector3 = Eigen::Matrix<Real, 1, 3>;
using Vector9r = Eigen::Matrix<Real, 9, 1, Eigen::DontAlign>;
using Matrix9r = Eigen::Matrix<Real, 9, 9, Eigen::DontAlign>;
using Matrix32r = Eigen::Matrix<Real, 3, 2, Eigen::DontAlign>;
using Matrix34r = Eigen::Matrix<Real, 3, 4, Eigen::DontAlign>;
using Matrix43r = Eigen::Matrix<Real, 4, 3, Eigen::DontAlign>;
using Matrix23r = Eigen::Matrix<Real, 2, 3, Eigen::DontAlign>;
using Matrix39r = Eigen::Matrix<Real, 3, 9, Eigen::DontAlign>;

using Matrix1212r = Eigen::Matrix<Real, 12, 12, Eigen::DontAlign>;
using Vector12r = Eigen::Matrix<Real, 12, 1, Eigen::DontAlign>;

template<typename T>
using MatSparse = Eigen::SparseMatrix<T>;

template<typename T, unsigned int N, unsigned int M>
using MatS = Eigen::Matrix<T, N, M, Eigen::DontAlign>;

template<typename T, unsigned int N>
using VecS = Eigen::Matrix<T, N, 1, Eigen::DontAlign>;

template <typename T>
using MatD = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using VecD = Eigen::Matrix<T, Eigen::Dynamic, 1>;

#endif //COMMON_H
