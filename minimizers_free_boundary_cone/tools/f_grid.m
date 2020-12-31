function fun = f_grid(beta,t,N,one,err)

if size(beta,2) > 1
    beta = beta.';
end
if size(t,1) > 1
   t = t.'; 
end

half = one/2;
u = 1-t;

a_old = -half*beta;

fun = 1+a_old*u;
for n = 2:N
    
    an = half*a_old.*(1-(n+beta)/n^2);
    fun = fun + an*u.^n;
    a_old = an;

end

fun = fun + iv(-err,err); 
    
