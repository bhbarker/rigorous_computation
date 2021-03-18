function x0 = newton(fun,x0,iter)
% =================================================================
%Purpose of the function:
% returns the root of the function with Newton's methond
% Parameters:
% fun: function 
% x0: initial guess
% iter: number of interation
% =================================================================
for j = 1:iter
    [f,f_der] = fun(x0);
   
    x0 = x0-f_der\f;
end
