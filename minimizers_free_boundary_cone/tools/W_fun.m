function [W,W_x,W_y] = W_fun(x,y,beta,mu,h_fun,g_fun,h_fun_der,g_fun_der,iv_fun)

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

one = iv_fun(1);

R = sqrt(x.^2+y.^2);
temp_h = h_fun(x./R,beta);
temp_g = g_fun(x./R,beta);

W = R.*temp_h+mu*temp_g./R.^(one/2);

x_r =x./R;
y2 = y.^2;
temp_h = h_fun(x_r,beta);
temp_g = g_fun(x_r,beta);

temp_h_der = h_fun_der(x_r,beta);
temp_g_der = g_fun_der(x_r,beta);

W_x = x_r.*temp_h+R.*temp_h_der.*(y2./R.^3) ...
-(mu/iv_fun(2))*x./(R.^(iv_fun(5)/2)).*temp_g+ ...
mu*temp_g_der.*(y2./R.^3)./R.^(one/2);

W_y = y.*temp_h./R+R.*temp_h_der.*(-x.*y)./R.^3- ...
(mu/iv_fun(2))*y.*temp_g./R.^(iv_fun(5)/2)+ ...
mu*temp_g_der.*(-x.*y)./R.^(iv_fun(7)/2);

















