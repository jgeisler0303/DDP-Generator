function out= optDefSimple

syms u1 x1 x2 x3 x4 dx real

A_= expm(diag(ones(4, 1), 1)*dx);
A= A_(1:4, 1:4);
B= A_(1:4, 5);

x_= [x1, x2, x3, x4]';
u_= [u1];

out.f= A*x_ + B*u_;
out.F= x1*x1;
out.L= x1*x1 + u1*u1;

out.x_= [x1 x2 x3 x4];
out.u_= [u1];
out.params= [dx];