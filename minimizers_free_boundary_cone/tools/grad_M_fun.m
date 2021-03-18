function out = grad_M_fun(x,y,beta,zk,cf_legendre,iv_fun)
% =================================================================
%Purpose of the function:
%calculate the gradient of M
% Parameters:
% x: independent variable
% y: independent variable
% beta: independent variable
% zk: coefficients of f
% cf_legendre: coefficients of legendre polynomial
% iv_fun: option to use double or interval arithmetic
% =================================================================

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

% one = iv_fun(1); % interval enclosure of 1
N = length(zk); % the number of coefficients 
k = iv_fun(0:1:N-1); % index


bk = (-1 + sqrt(1+(1./beta)*(8*(k+1).*k)))/2; % The powers of R

% The Legendre polynomials
M = iv_fun(zeros(length(x),1));


R = sqrt(x.^2+y.^2);

dom = x./R;

sm1 = iv(0);
sm2 = iv(0);
for k = 0:N-1
    lk = iv_fun(zeros(size(R)));
    lk_der = iv_fun(zeros(size(R)));
    for n = 0:k
        lk = lk + cf_legendre(k+1,n+1)*dom.^n;
        lk_der = lk_der + n*cf_legendre(k+1,n+1)*dom.^(n-1);
    end

    M = M + (zk(k+1)*R.^bk(:,k+1).*lk);
    
    temp1 = zk(k+1)*bk(:,k+1).*R.^(bk(:,k+1)-1).*lk;
    temp2 = zk(k+1)*R.^bk(:,k+1).*lk_der;
    
    sm1 = sm1 + temp1;
    sm2 = sm2 + temp2;


    
    
end

out = (beta/2).*sm1.^2+(1-(x./R).^2).*(1./R.^2).*sm2.^2;










