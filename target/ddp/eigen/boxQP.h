#ifndef BOXQP_H
#define BOXQP_H

#include "iLQG.hpp"
#include <Eigen/Dense>

using namespace Eigen;

static const int nQP= N_U;

// https://eigen.tuxfamily.org/dox/classEigen_1_1Ref.html
int boxQP(const Ref<const MatrixUU> &H, const Ref<const VectorU> &g, const Ref<const VectorU> &lower, const Ref<const VectorU> &upper, Ref<VectorU> x, Ref<VectorU_int> is_clamped, int &n_free, LLT<MatrixUU_dyn, Upper> &llt);
// int boxQP(MatrixUU &H, const Ref<const VectorU> &g, const Ref<const VectorU> &lower, const Ref<const VectorU> &upper, Ref<VectorU> x, Ref<VectorU_int> is_clamped, int &n_free, Ref<MatrixUU_dyn> lltHfree);

#endif /* BOXQP_H */
