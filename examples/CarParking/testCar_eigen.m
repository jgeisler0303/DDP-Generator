randn ("state", 0);
%% parameters
p.d= 2;
p.h= 0.03;

p.pf= [.01 .01 .01  1];
p.cf= [ .1  .1   1  .3];
p.cu= 1e-2*[1 .01];
p.cx= 1e-3*[1  1];
p.px= [.1 .1];
p.limW= [-.5 .5];
p.limA= [-2  2];

%% initial conditions
T= 500;              % horizon
t= (1:T+1)*p.h;
x0= [1;1;pi*3/2;0];   % initial state
u0= .1*randn(2,T);    % initial controls

Op.max_iter= 200;

tic
[success_, x_, u_, cost_, cx_, cu_, cxx_, cuu_, cxu_, fx_, fu_, fxx_, fuu_, fxu_]= ddpCar_eigen3(x0, u0, p, Op);
toc

plotOptCar(t, x_, u_, [], [], p)
