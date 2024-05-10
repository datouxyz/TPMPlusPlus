#pragma once
#include <Eigen/Dense>
#include "ParamaterDef.h"

MAT_TYPE RDivide(MAT_TYPE &A, MAT_TYPE &B);

MAT_TYPE_F RDivide(MAT_TYPE_F &A, MAT_TYPE_F &B);

MAT_TYPE_F LDivide(MAT_TYPE_F &A, MAT_TYPE_F &B);

MAT_TYPE LDivide(MAT_TYPE &A, MAT_TYPE &B);
