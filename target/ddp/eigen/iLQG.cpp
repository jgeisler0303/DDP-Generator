// C implementation of iLQG algorithm from http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization by Yuval Tassa
// Copyright (c) 2016 Jens Geisler
//
// BIBTeX:
// @INPROCEEDINGS{
// author={Tassa, Y. and Mansard, N. and Todorov, E.},
// booktitle={Robotics and Automation (ICRA), 2014 IEEE International Conference on},
// title={Control-Limited Differential Dynamic Programming},
// year={2014}, month={May}, doi={10.1109/ICRA.2014.6907001}}

#include "iLQG.hpp"


tOptSet o= optSet();
    