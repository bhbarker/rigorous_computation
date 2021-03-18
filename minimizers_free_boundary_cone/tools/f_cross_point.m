function [f,J] = f_cross_point(x,y,beta,zk_standard,cf_legendre_standard_basis,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun)
%  Purpose:
%  Finding the croos point of function M and W

% function M and partial derivatives
[M,M_x,M_y] = M_fun(x,y,beta,zk_standard,cf_legendre_standard_basis,iv_fun); 
% function W and partial derivatives
[W,W_x,W_y] = W_fun(x,y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun);

f = mid([M;W]);
% Jacobian matrix
J = mid([M_x, M_y; W_x, W_y]);

