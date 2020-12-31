function out = h_grid(beta,t,f_fun,f_der_fun,x0_fun,err)

if size(beta,2) > 1
    beta = beta.';
end
if size(t,1) > 1
   t = t.'; 
end

temp = f_der_fun(x0_fun(beta),beta).*sqrt(1-x0_fun(beta).^2);

out = f_fun(beta,t)./repmat(temp,1,length(t));

out = out + iv(-err,err);