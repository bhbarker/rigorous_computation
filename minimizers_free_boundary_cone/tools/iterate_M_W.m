function [X,Y] = iterate_M_W(x,y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun,zk_standard,cf_legendre_standard_basis)



[M,M_x,M_y] = M_fun(x,y,beta,zk_standard,cf_legendre_standard_basis,iv_fun); 
[W,W_x,W_y] = W_fun(x,y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun);

[M0,M0_x,M0_y] = M_fun(mid(x),mid(y),beta,zk_standard,cf_legendre_standard_basis,iv_fun); 
[W0,W0_x,W0_y] = W_fun(mid(x),mid(y),beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun);

dt = 1./(M_x.*W_y-W_x.*M_y);

ind = find(sup(dt).*inf(dt)<0);
if ~(isempty(ind))
   error('Determinant includes zero.'); 
end

X = intersect(x,mid(x) - (W_y.*M0-M_y.*W0).*dt);
Y = intersect(y,mid(y) - (M_x.*W0-W_x.*M0).*dt);


