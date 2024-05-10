#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

enum BsxOperator
{
	BsxOperator_times,
	BsxOperator_plus,
	BsxOperator_minus,
	BsxOperator_divide,
	BsxOperator_rounded_divide,
	BsxOperator_MAX = 0xffffffff
};

#define MAT_TYPE Eigen::MatrixXd
#define MAT_TYPE_F Eigen::MatrixXf

#define MAT_SPARSE_TYPE Eigen::SparseMatrix<double>

#define VECTOR_TYPE Eigen::VectorXd

#define VECTOR_TYPE_F Eigen::VectorXf

#define TRIPLET_TYPE Eigen::Triplet<double>
typedef std::vector<TRIPLET_TYPE> TRIPLETS;

typedef std::vector<VECTOR_TYPE> VECTOR_TYPES;

typedef std::vector<VECTOR_TYPES> VECTOR_TYPES_VEC;


enum RegularMode
{
	RegularMode_BentEnergy,
	RegularMode_MembraneEnergy,
	RegularMode_MAX = 0xffffffff
};

