function out = f(t,beta,N,one,err)
% =================================================================
%Purpose of the function:
%evaluate function f
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
if sty > stx % make beta a column vector
   beta = beta.';
end

%Power series of f round (1-t)
half = one/2;
u = 1-t;

a_old = -half*beta;

out = 1+u.*a_old;
for n = 2:N
    
    an = half*a_old.*(1-(n+beta)/n^2);
    out = out + u.^n.*an;
    a_old = an;

end



if nargin == 5
    
    if err == 0
        return
    end
    out = out + iv(-err,err); 
    
end
