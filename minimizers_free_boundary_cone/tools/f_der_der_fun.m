function out = f_der_der_fun(t,beta,N,one,err)
% =================================================================
%Purpose of the function:
%evaluate the second derivative of function f
% Parameters:
% t: independent variable
% beta: independent variable
% N: number of terms for the power series
% one: option to use double arithmetic and interval arithmetic
% err: values to evaluate
% =================================================================
[stx,sty] = size(t);
if sty > stx % make t a column vector
   t = t.';
end

[stx,sty] = size(beta);
if stx > sty % make beta a row vector
   beta = beta.';
end

half = one/2;
u = 1-t;

out = 0*one;
a_old = -half*beta;
for n = 2:N
    
    an = half*a_old.*(1-(n+beta)/n^2);
    out = out + n*(n-1)*u.^(n-2)*an;
    a_old = an;

end

if nargin == 5
   out = out + iv(-err,err); 
end


