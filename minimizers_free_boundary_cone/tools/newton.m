function x0 = newton(fun,x0,iter)

for j = 1:iter
    [f,f_der] = fun(x0);
   
    x0 = x0-f_der\f;
end
