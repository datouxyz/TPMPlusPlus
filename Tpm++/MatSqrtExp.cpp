#include <Eigen/Dense>
#include <Eigen/MatrixFunctions>
#include "ParamaterDef.h"

MAT_TYPE Sqrtm(MAT_TYPE m)
{
	return m.sqrt();
}
VECTOR_TYPE Sqrtm(VECTOR_TYPE m)
{
	return m.sqrt();
}

MAT_TYPE Expm(MAT_TYPE m)
{
	return m.exp();
}
MAT_TYPE Logm(MAT_TYPE m)
{
	return m.log();
}