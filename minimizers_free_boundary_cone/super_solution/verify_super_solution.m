function verify_super_solution(file_name,c_int,small_interval,zk_data,zf)
%{
Takes as input:
 
file_name (string) - Data is saved as "verified_" concatenated with file_name.
c_int (2x1 array of intervals) - Contains two intervals, the left and right
end point of the interval of c-values that will be verified.
small_interval (interval) - The cross point will be rigorously verified by
initializing interval Newton's method with an approximation of the cross
point plus the small_interval.
zk_data (array of intervals) - The coefficients for M.
zf (structure) - The structure that contains data about the zeros of f.
%}

t_start = tic;

% Make sure all are indeed intervals
zk_standard = iv(zk_data.zk_standard); % make coefficients into intervals
iv_fun = @(x)iv(x); % do rigorous computations
mu = iv(zk_data.mu); % rigorous mu

% beta intervals to check
beta_temp = 2./(1+c_int.^2);
beta_a = min(inf(beta_temp));
beta_b = max(sup(beta_temp));
beta_lin = linspace(beta_a,beta_b,1001);
beta = iv(beta_lin(1:end-1),beta_lin(2:end)); 


%{
    Get interpolation error bounds for f(t,beta) and h(t,beta)
%}

% Interval of which t is an element.
a_t = iv('-0.7'); % left end point of t interval
b_t = iv(1); % right end point of t interval

% rho giving stadium standard_basis_from_chebyshev
rho_t = iv('1.5');
rho_beta = iv('1.5');

[h0_fun_der,h0_fun,h0_grid] = f_and_h(zf,iv_fun,a_t,b_t,rho_beta,rho_t);

[g0_fun,g0_fun_der,g0_grid] = interp_g(iv('-0.25'),iv('2'),a_t,b_t,rho_beta,rho_t);


cf_legendre_standard_basis = legendre_polys(length(zk_standard),'on');

t = linspace(-1,1,10001);
t = iv(t(1:end-1),t(2:end));
N = length(zk_standard); % the number of coefficients 
M = iv(zeros(1,length(t)));
for k = 0:N-1
    lk = iv(zeros(1,length(t)));
    for n = 0:k
        lk = lk + cf_legendre_standard_basis(k+1,n+1)*t.^n;
    end
    M = M + zk_standard(k+1)*lk; 
end
ind = find(sup(M)>0);
tc = t(ind);

verify_condition_for_super_solution(tc,beta,zk_standard, ...
    cf_legendre_standard_basis,mu,g0_grid,h0_grid);

MW = @(x,y) iterate_M_W(x,y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun,zk_standard,cf_legendre_standard_basis);

X = iv(zk_data.x0(1)*ones(length(beta),1));
Y = iv(zk_data.x0(2)*ones(length(beta),1));
for k = 1:10
    
    [M0,M0_x,M0_y] = M_fun(X,Y,beta,zk_standard,cf_legendre_standard_basis,iv_fun); 
    [W0,W0_x,W0_y] = W_fun(X,Y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun);
    dt = 1./(M0_x.*W0_y-W0_x.*M0_y);
    X = mid(X - (W0_y.*M0-M0_y.*W0).*dt);
    Y = mid(Y - (M0_x.*W0-W0_x.*M0).*dt);    
    
    
%     M0(1)
%     W0(1)
%     fprintf('\n\n');

end


X = X+small_interval;
Y = Y+small_interval;
for j = 1:10
    
    [X,Y] = MW(X,Y);
end


X(1)
Y(1)
%                             

% [M,M_x,M_y] = M_fun(X,Y,beta,zk_standard,cf_legendre_standard_basis,iv_fun); 
% [W,W_x,W_y] = W_fun(X,Y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der,iv_fun);

cond1 = grad_W_fun(X,Y,beta,mu,h0_fun,g0_fun,h0_fun_der,g0_fun_der);
cond2 = grad_M_fun(X,Y,beta,zk_standard,cf_legendre_standard_basis,iv_fun);

max_condition_1 = max(sup(cond1));
max_condition_2 = max(sup(cond2));

if sum(isnan(max_condition_1)) > 0
   error('Failed to verify'); 
end
if sum(isnan(max_condition_2)) > 0
   error('Failed to verify'); 
end

if max_condition_1 >= 1
    error('Failed to verify supersolution.');
end
if max_condition_2 >= 1
    error('Failed to verify supersolution.');
end
 
R = bound_R_for_M_zero(beta,zk_standard);

verify_unique_cross_point(R,Y,beta,zk_standard,cf_legendre_standard_basis);

t_stop = toc(t_start);

data.zk_data = zk_data;
data.c_int = c_int;
data.verified = 'yes';
data.small_interval = small_interval;
data.R = R;
data.run_time = t_stop;

curr_dir = cd;
cd('../data');
save(['verified_',file_name],'data');
cd(curr_dir);







