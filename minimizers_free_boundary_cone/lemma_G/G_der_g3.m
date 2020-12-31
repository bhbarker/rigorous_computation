function out = G_der_g3(t,beta,f0_fun,f1_fun,f2_fun,g0_fun,g1_fun,g2_fun)

f0 = f0_fun(t,beta);
f1 = f1_fun(t,beta);
f2 = f2_fun(t,beta);
g0 = g0_fun(t,beta);
g1 = g1_fun(t,beta);
g2 = g2_fun(t,beta);


[stx,sty] = size(beta);
if stx > sty
   beta = beta.';
end

beta_mat = repmat(beta,size(f0,1),1);

[stx,sty] = size(t);
if sty > stx
   t = t.';
end
t_mat = repmat(t,1,size(f0,2));




con = iv(9)/8;

out = 2*con*beta_mat.*f0.*f1.*g0.^3 - ...
    2*t_mat.*(f1.*g0-f0.*g1).^2.*g0+ ...
    2*(1-t_mat.^2).*(f1.*g0-f0.*g1).* ...
    (f2.*g0.^2-f1.*g1.*g0 - f0.*g2.*g0+f0.*g1.^2); 






