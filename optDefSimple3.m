function out= optDefSimple3

syms u1 u2 x1 x2 x3 x4 dx e_ ref_k real

A_= expm(diag(ones(4, 1), 1)*dx);
A= A_(1:4, 1:4);
B= A_(1:4, 5);

x_= [x1, x2, x3, x4]';
u_= [u1];

out.f= A*x_ + B*u_ + [1 0 0 0]'*u2;
out.F= (x1-ref_k)*(x1-ref_k);
out.L= (x1-ref_k)*(x1-ref_k) + sqrt(u1*u1+e_) + exp(u2*u2);

out.x_= [x1 x2 x3 x4];
out.u_= [u1 u2];
out.params= [dx e_ ref_k];