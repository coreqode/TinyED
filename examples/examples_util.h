#ifndef EXAMPLES_UTIL_H
#define EXAMPLES_UTIL_H

#include <set>

#include "common.h"
#include "filesystem"
#include <fstream>

// Taken from trimesh library
namespace trimesh
{
    typedef long index_t;

    struct edge_t
    {
        index_t v[2];

        index_t& start() { return v[0]; }
        const index_t& start() const { return v[0]; }

        index_t& end() { return v[1]; }
        const index_t& end() const { return v[1]; }

        edge_t()
        {
            v[0] = v[1] = -1;
        }

        edge_t(index_t i, index_t j)
        {
            v[0] = i;
            v[1] = j;
        }

        bool operator<(const edge_t& other) const {
            if (v[0] != other.v[1])
                return v[0] < other.v[0];
            return v[1] < other.v[1];
        }
    };

    struct triangle_t
    {
        index_t v[3];

        index_t& i() { return v[0]; }
        const index_t& i() const { return v[0]; }

        index_t& j() { return v[1]; }
        const index_t& j() const { return v[1]; }

        index_t& k() { return v[2]; }
        const index_t& k() const { return v[2]; }

        triangle_t()
        {
            v[0] = v[1] = v[2] = -1;
        }

        triangle_t(index_t i, index_t j, index_t k)
        {
            v[0] = i;
            v[1] = j;
            v[2] = k;
        }
    };
}

using namespace TinyED;

struct Element
{
    std::variant<trimesh::edge_t, trimesh::triangle_t> elementData;
};

std::tuple<std::vector<trimesh::triangle_t>, Mat3N, Mat2N> generateTesselatedSquare(Vec3 center, size_t n, Real sideLength)
{
    std::vector<trimesh::triangle_t> triangles;
    triangles.reserve(n * n * 2);
    trimesh::index_t verticesPerRow = n + 1;
    Mat3N positions = center.replicate(1, verticesPerRow * verticesPerRow);
    Mat2N uvCoords = Mat2N::Zero(2, verticesPerRow * verticesPerRow);
    Real tileLength = sideLength / n;
    Real halfLength = sideLength / 2;
    for (trimesh::index_t i = 0; i < n; i++)
    {
        for (trimesh::index_t j = 0; j < n; j++)
        {
            trimesh::index_t prev = i * verticesPerRow;
            triangles.push_back({
                j + prev,
                j + prev + verticesPerRow,
                j + prev + verticesPerRow + 1,
            });
            triangles.push_back({
                j + prev + 1,
                j + prev,
                j + prev + verticesPerRow + 1,
            });
        }
    }
    for (size_t i = 0; i <= n; i++)
    {
        for (size_t j = 0; j <= n; j++)
        {
            uvCoords.col(i * verticesPerRow + j) = Vec2(i * sideLength, j * sideLength) / n;
            positions.col(i * verticesPerRow + j) += Vec3({ j * tileLength - halfLength, 0, i * tileLength - halfLength });
            // positions.col(i * verticesPerRow + j) += Vec3({ 0, j * tileLength - halfLength, i * tileLength - halfLength });
            // positions.col(i * verticesPerRow + j) += Vec3({ i * tileLength - halfLength, j * tileLength - halfLength, 0 });
        }
    }
    return { triangles, positions, uvCoords };
}

// Computes the rest shapes for triangles
std::vector<Mat2> computeRestShape(const Mat3N& positions, const Mat2N& UVs, const std::vector<trimesh::triangle_t>&
triangles)
{
    const trimesh::index_t triangleCount = triangles.size();
    std::vector<Mat2> restMatrices(triangleCount);

    for (trimesh::index_t t = 0; t < triangleCount; t++)
    {
        const unsigned int v0 = triangles[t].v[0];
        const unsigned int v1 = triangles[t].v[1];
        const unsigned int v2 = triangles[t].v[2];


        Mat2 UVDiffs; UVDiffs.setZero();
        UVDiffs.col(0) = UVs.col(v1) - UVs.col(v0);
        UVDiffs.col(1) = UVs.col(v2) - UVs.col(v0);

        restMatrices[t] = UVDiffs.inverse() ;
    }
    return restMatrices;
}

std::vector<Real> computeRestVolume(const Mat2N& UVs, const std::vector<trimesh::triangle_t>& triangles)
{
    const trimesh::index_t triangleCount = triangles.size();
    std::vector<Real> restVolumes(triangleCount);

    for (trimesh::index_t t = 0; t < triangleCount; t++)
    {
        const unsigned int v0 = triangles[t].v[0];
        const unsigned int v1 = triangles[t].v[1];
        const unsigned int v2 = triangles[t].v[2];

        Vec2 uv0 = UVs.col(v0);
        Vec2 uv1 = UVs.col(v1);
        Vec2 uv2 = UVs.col(v2);

        Mat2 uv;
        uv.col(0) = uv1 - uv0;
        uv.col(1) = uv2 - uv0;

        restVolumes[t] = 0.5 * std::abs(uv.determinant());
    }
    return restVolumes;
}

// Computes the rest shapes for edges
std::vector<Real> computeRestShape(const Mat3N& positions, const Mat2N& UVs, const std::vector<trimesh::edge_t>& edges)
{
    std::vector<Real> restLengths;

    return restLengths;
}


VecN computeVertexMasses(const Mat3N& positions, const std::vector<trimesh::triangle_t>& triangles, const Real density)
{
    const trimesh::index_t vertexCount = positions.cols();
    VecN masses(vertexCount);
    masses.setZero();

    for (trimesh::triangle_t t : triangles)
    {
        const Real area = (positions.col(t.j()) - positions.col(t.i())).cross(positions.col(t.k()) - positions.col(t.i())
        ).norm() / 2;
        const Real massContrib = area * density;

        masses(t.i()) += massContrib / 3;
        masses(t.j()) += massContrib / 3;
        masses(t.k()) += massContrib / 3;
    }

    return masses;
}

template <unsigned int N>
void computeMaterialForce(const VecS<3*N> &x,
                          const Material &mat,
                          VecS<3*N> &forceLocal)
{
    forceLocal.setZero();

    switch (mat.type)
    {
    case MaterialType::SPRING:
            if constexpr (N == 2)
                SpringModel::getForce(x, std::get<SpringData>(mat.data), forceLocal);
            break;
    case MaterialType::STVK:
            if constexpr (N == 3)
                StVKModel::getForce(x, std::get<FEMTriData>(mat.data), forceLocal);
            break;
    case MaterialType::NEOHOOKEAN:
            throw std::runtime_error( "Error: NEOHOOKEAN material not yet implemented.");
    default:
            throw std::runtime_error("Error: Unknown material type.");
    }
}


template <unsigned int N>
void computeMaterialJacobian(const VecS<3 * N>& x, const Material& mat, MatS<3 * N, 3 * N>& jacobianLocal)
{
    jacobianLocal.setZero();

    switch (mat.type)
    {
    case MaterialType::SPRING:
            if constexpr (N == 2)
                SpringModel::getForceJacobian(x, std::get<SpringData>(mat.data), jacobianLocal);
            break;
    case MaterialType::STVK:
            if constexpr (N == 3)
                StVKModel::getForceJacobian(x, std::get<FEMTriData>(mat.data), jacobianLocal);
            break;
    case MaterialType::NEOHOOKEAN:
            throw std::runtime_error(
                "Error: NEOHOOKEAN material not implemented.");
    default:
            throw std::runtime_error( "Error: Unknown or unimplemented material type.");
    }
}



void computeAllForces(const Mat3N& positions, const Mat3N& velocities, const std::vector<Element>&
elements, const std::vector<Material>& matData, Mat3N& forces)
{
    forces.setZero();

    for (size_t i = 0; i < elements.size(); i++)
    {
        std::visit([&](auto&& el) {
            using T = std::decay_t<decltype(el)>;
            constexpr size_t N = std::is_same_v<T, trimesh::triangle_t> ? 3 : 2;

            VecS<3 * N> x;
            x.setZero();

            if constexpr (N == 2)
            {
                x.segment(0, 3) = positions.col(el.v[0]);
                x.segment(3, 3) = positions.col(el.v[1]);
            }
            else if constexpr (N == 3)
            {
                x.segment(0, 3) = positions.col(el.v[0]);
                x.segment(3, 3) = positions.col(el.v[1]);
                x.segment(6, 3) = positions.col(el.v[2]);
            }

            VecS<3 * N> forceLocal;
            forceLocal.setZero();
            computeMaterialForce<N>(x, matData[i], forceLocal);

            if constexpr (N == 2)
            {
                forces.col(el.v[0]) += forceLocal.segment(0, 3);
                forces.col(el.v[1]) += forceLocal.segment(3, 3);
            }
            else if constexpr (N == 3)
            {
                forces.col(el.v[0]) += forceLocal.segment(0, 3);
                forces.col(el.v[1]) += forceLocal.segment(3, 3);
                forces.col(el.v[2]) += forceLocal.segment(6, 3);
            }
        }, elements[i].elementData);
    }
}

void computeAllForceJacobians(const Mat3N& positions, const Mat3N& velocities, const std::vector<Element>&
                              elements, const std::vector<Material>& matData, MatN& forceJacobians)
{
    forceJacobians.setZero();
    for (size_t i =0; i < elements.size(); i++)
    {

        std::visit([&](auto&& el)
        {
            using T = std::decay_t<decltype(el)>;
            constexpr size_t N = std::is_same_v<T, trimesh::triangle_t> ? 3 : 2;

            VecS<3 * N> x;
            x.setZero();

            if constexpr (N == 2)
            {
                x.segment(0, 3) = positions.col(el.v[0]);
                x.segment(3, 3) = positions.col(el.v[1]);
            }
            else if constexpr (N == 3)
            {
                x.segment(0, 3) = positions.col(el.v[0]);
                x.segment(3, 3) = positions.col(el.v[1]);
                x.segment(6, 3) = positions.col(el.v[2]);
            }

            Mat<3*N, 3*N> forceJacobianLocal;
            forceJacobianLocal.setZero();
            computeMaterialJacobian<N>(x, matData[i], forceJacobianLocal);

            for (size_t j = 0; j < N; j++)
            {
                for (size_t k =0; k < N; k++)
                {
                    forceJacobians.template block<3,3>(el.v[j] * 3, el.v[k] * 3) += forceJacobianLocal
                    .template block<3,3>
                    (j * 3,
                     k * 3);
                }
            }
        }, elements[i].elementData);
    }
}


void saveSheet(const std::string fileName, const Mat3N& positions, const std::vector<trimesh::triangle_t>& triangles)
{
    if (!std::filesystem::exists("./output"))
    {
        std::filesystem::create_directory("./output");
    }
    std::ofstream file("./output/" + fileName + ".obj");
    if (file.is_open())
    {
        for (size_t i = 0; i < positions.cols(); i++)
        {
            file << "v " << positions(0, i) << " " << positions(1, i) << " " << positions(2, i) << std::endl;
        }
        for (size_t i = 0; i < triangles.size(); i++)
        {
            file << "f " << triangles[i].i() + 1 << " " << triangles[i].j() + 1 << " " << triangles[i].k() + 1 << std::endl;
        }
        file.close();
    }
}

void advanceStep(Mat3N& positions, Mat3N& velocities, const Real timeStep, const std::vector<Element>&
elements, const std::vector<Material>& matData, const VecN& invMlumped, const VecN& masses, const Vec3& gravity)
{
        const Mat3N vPrev = velocities;
        const Mat3N xPrev = positions;

        constexpr size_t m_Iter = 2;
        constexpr Real m_Eps = 1e-6;

        const unsigned int N = positions.cols();

        Mat3N forces(3, N); forces.setZero();
        MatN dfDx(N*3, N*3); dfDx.setZero();
        MatN kron(N*3,N*3); kron.setIdentity();

        for (size_t iter = 0; iter < m_Iter; iter++)
        {
            positions = xPrev + timeStep * velocities;

            forces.setZero();
            computeAllForces(positions, velocities, elements, matData, forces);

            // add external gravity force
            for (size_t i = 0; i < forces.cols(); i++)
            {
                forces.col(i) += masses(i) * gravity;
            }

            const VecN residual = FlattenMatrixN3(velocities) - FlattenMatrixN3(vPrev) - timeStep * (invMlumped.asDiagonal()
            * FlattenMatrixN3(forces));

            if (residual.norm() < m_Eps)
                break;

            dfDx.setZero();
            computeAllForceJacobians(positions, velocities, elements, matData, dfDx);

            const MatN H = kron - timeStep * timeStep * (invMlumped.asDiagonal() * dfDx);
            VecN delta_v; delta_v.setZero();

            delta_v = H.partialPivLu().solve(-residual);
            velocities += GroupMatrixN3(delta_v);
        }
        positions = xPrev + timeStep * velocities;
}


#endif //EXAMPLES_UTIL
