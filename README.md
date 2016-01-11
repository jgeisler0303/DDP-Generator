Introduction
============
This is a highly efficient C implementation of [Yuval Tassa's iLQG algorithm](http://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization) (see also his [paper](https://homes.cs.washington.edu/~todorov/papers/TassaICRA14.pdf)).

Efficiency is achieved using the following features:
* Problem-specific code (system function, cost function and their derivatives) is generated from an analytic problem description in the [Maxima language](http://maxima.sourceforge.net/)
* Auxillari variables and their derivatives used in the system function and cost function are reused avoiding unneccesary recalculations
* For symmetric matrices only the upper triangle is calculated

For the car-parking example on my computer this yields a speed improvement of 8ms per iteration vs. 1500ms for the original MATLAB implementation.

Everything is still very raw, untested and undocumented. But if you are interested I will be more than happy to supply addition information and help with problems. In this case, please just open new a github issue in my repo.

Code Generation
===============
The problem-specific code is generated from a problem definition in the [Maxima language](http://maxima.sourceforge.net/) using the maxima package `gentran`. Unfortunately, the gentran package as supplied by the current maxima release is not fully functional. Therefore you will have to download the package from my [my repo](https://github.com/jgeisler0303/maxima) and replace the files in the `maxima/share/contrib/gentran` folder of your maxima installation with the files from mym repo. On my Linux computer I got gentran only to work using [Steel Bank Common Lisp (SBCL)](http://www.sbcl.org/) implementation. Unfortunately, this meant I had to recompile maxima from source. On Windows, SBCL is the default Lisp implementation supplied with the [wxMaxima](http://andrejv.github.io/wxmaxima/) installation.

Problem Description
-------------------
Your optimal control problem has to be defined in a maxima batch file, i.e. a text file with `.mac` suffix containing maxima expressions. In this file at least six variables have to be defined:
* `x`: a list of symbols for the state names (e.g. `x: [state1, state2, state3]`)
* `u`: a list of symbols for the system inputs (e.g. `x: [in1]`)
* `f`: the time discrete state transition function as an array index by the state name symbols. E.g.
```
    f[state1]: state2*ts;
    f[state2]: state3*ts+offset_parameter;
    f[state3]: in1*ts;
```
* `L`: the running cost function
* `F`: the final cost function (may not depend on inputs)
* `h`: an array of input contraint functions. For every `h[i]` it is assumed that `h[i]<0` and every constraint is enforced in order of ascending index of h. Each `h[i]` may only depend on one input with positive or negative one as coefficient.

The symbols for states and inputs my not be defined. In the expressions of `f`, `L`, `F` and `h` any undefined symbol appart from state and input symbols is considered to be a parameter. Any defined symbol is considered to be an auxillari value. Auxillari values and their derivatives are calculated separately and before their use in `f`, `L`, `F` or `h` or their derivatives. Auxillari values my depend on other auxillari values, states, inputs or parameters. If you use auxillaries, make sure to prefix them with an apostrophe. Otherwise, the auxillari definition will be substitued in to your expression and you end up without the benefits of auxillaries.

Getting the Example to Run
==========================
As an example the definition of the car parking problem from Yuval Tassa's iLQG implementation can be found in the examples folder. This definition can be compiled into an [Octave Mex function](https://www.gnu.org/software/octave/doc/interpreter/Getting-Started-with-Mex_002dFiles.html) using the `make_iLQG.sh` shell script. Assuming you have Maxima, Octave and the Octave `mkoctfile` package installed, simply go to the `CarParking` folder and type `../../make_iLQG.sh optDefCar`. Upon success, a subfolder named `Car_gen_files` is created, containing the generated problem specific code and an empty `err.log` file. In the `CarParking` foleder, the compiled `iLQGCar.mex` file is created. Now you can start Octave and run the 'testCar` demo script.

For Windows users there is yet no make script. But, if you are interested, I could write a script that can be run in Octave or MATLAB.

The script currently compiles a Mex-function for Octave. But since the Octave API is an almost perfect clone of the MATLAB API, a port to MATLAB only requires an additional `#include "matrix.h` and compilation using the MATLAB `mex`script instead of `mkoctfile`.

The expamle will not always converge to the same result as the original by Yuval Tassa. This is mostly due to the fact that the example is very sensitive to small numerical differences. Most of these differences are introduced by the difference in the calculation of the derivatives: the generated code calculates analytic derivatives, while Yuval's example uses finite differences approximations. But in my investigations, already formulating the problem with and without auxillaries introduced enough differences to result in different convergence behaviour.

