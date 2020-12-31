function W_y = W_temp(x,y,beta,mu,h_fun,g_fun,h_fun_der,g_fun_der,iv_fun)

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

R = sqrt(x.^2+y.^2);


x_r =x./R;
temp_h = h_fun(x_r,beta);
temp_g = g_fun(x_r,beta);

temp_h_der = h_fun_der(x_r,beta);
temp_g_der = g_fun_der(x_r,beta);


W_y = R.^(iv(5)/2).*temp_h-x.*R.^(iv(3)/2).*temp_h_der-(mu/2)*R.*temp_g-x.*temp_g_der;
















