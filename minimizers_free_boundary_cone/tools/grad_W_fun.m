function out = grad_W_fun(x,y,beta,mu,h_fun,g_fun,h_fun_der,g_fun_der)

% make sure that x is a column vector
if size(x,2) > 1
   x = x.'; 
end

% make sure that y is a column vector
if size(y,2) > 1
   y = y.'; 
end

% make sure that beta is a column vector
if size(beta,2) > 1
   beta = beta.'; 
end

% one = iv_fun(1);

R = sqrt(x.^2+y.^2);
% temp_h = h_fun(x./R,beta);
% temp_g = g_fun(x./R,beta);

% W = R.*temp_h+mu*temp_g./R.^(one/2);

x_r =x./R;
% y2 = y.^2;
temp_h = h_fun(x_r,beta);
temp_g = g_fun(x_r,beta);

temp_h_der = h_fun_der(x_r,beta);
temp_g_der = g_fun_der(x_r,beta);


out = (beta/2).*(temp_h-(mu./(2*R.^(iv(3)/2))).*temp_g).^2 + ...
    (1-x_r.^2).*(1./R.^2).*(R.*temp_h_der+(mu./sqrt(R)).*temp_g_der).^2;
















