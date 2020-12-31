function out = f_der_der_fun(t,beta,N,one,err)

[stx,sty] = size(t);
if sty > stx
   t = t.';
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


