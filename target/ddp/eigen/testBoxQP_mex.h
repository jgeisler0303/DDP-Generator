#ifndef TESTBOXQP_MEX_H
#define TESTBOXQP_MEX_H

typedef Eigen::Matrix<double, N_U, N_U> MatrixUU;
typedef Eigen::Matrix<double, Dynamic, Dynamic, ColMajor, N_U, N_U> MatrixUU_dyn;
typedef Eigen::Matrix<double, N_U, 1> VectorU;
typedef Eigen::Matrix<double, Dynamic, 1, ColMajor, N_U, 1> VectorU_dyn;
typedef Eigen::Matrix<int32_t, N_U, 1> VectorU_int;


#endif /* TESTBOXQP_MEX_H */
