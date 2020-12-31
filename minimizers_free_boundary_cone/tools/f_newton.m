function [fun,fun_der] = f_newton(t,beta,N,one)


half = one/2;
u = 1-t;


a_old = -half*beta;

fun = 1+a_old.*u;
fun_der = a_old;
for n = 2:N
    
    an = half*a_old.*(1-(n+beta)/n^2);
    fun = fun + an.*u.^n;
    fun_der = fun_der + n*an.*u.^(n-1);
    a_old = an;

end


fun_der = -fun_der;
