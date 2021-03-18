function verify_condition_for_super_solution(t,beta,zk, ...
    cf_legendre,mu,g0_grid,h0_grid)

% =================================================================
%Purpose of the function:
%verify the inequality given in Lemma 5.1 (1)
% Parameters:
% t: independant variable
% beta: independant variable
% zk: coefficients of f
% mu: constant variable
% cf_legendre: function that produce legendre basis
% g0_grid: function g
% h0_grid: function h
% =================================================================


% make sure that beta is a column vector
if size(beta,2) > 1
   beta = beta.'; 
end

% make sure that t is a row vector
if size(t,1) > 1
    t = t.';
end

temp_h = h0_grid(beta,t);
temp_g = g0_grid(beta,t);

W = temp_h+mu*temp_g;

% Compute M
N = length(zk); % the number of coefficients 
M = iv(zeros(1,length(t)));
for k = 0:N-1
    lk = iv(zeros(1,length(t)));
    for n = 0:k
        lk = lk + cf_legendre(k+1,n+1)*t.^n;
    end
    M = M + zk(k+1)*lk; 
end

check = repmat(M,length(beta),1)-W;


if min(min(inf(check)))<0
    error('Failed to verify condition');
end

