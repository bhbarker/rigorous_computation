function [h0_fun_der,h0_fun,h0_grid] = f_and_h(zf,iv_fun,a_t,b_t,rho_beta,rho_t)




cf_standard = standard_basis_from_chebyshev(zf.cf_p,iv_fun);

x0_fun = @(beta) eval_standard(cf_standard,zf.a,zf.b,beta) ...
    +iv(-zf.root_error_bound,zf.root_error_bound);


% (rho*exp(1i*theta)+exp(-1i*theta)/rho)/2
theta = linspace(0,sup(iv('pi')*2),101);
theta = iv(theta(1:end-1),theta(2:end));
E_beta_tilde = (rho_beta*exp(1i*theta)+1./(rho_beta*exp(1i*theta)))/2;
E_beta = (zf.a+zf.b)/2+(zf.b-zf.a)*E_beta_tilde/2;

% The square of the root of f given beta
y = x0_fun(E_beta).^2;


% Find the max real part of x0^2(beta) to show sqrt(1-x0^2(beta)) is
% analytic.
max_real_x0_squared = max(sup(real(y)));

if max_real_x0_squared >= 1
   error('Proof failed because sqrt{1-x_0^2(beta)} is not analytic on and inside the ellipse.')
end

% Get interpolation bounds on f' and on sqrt(1-x0^2(beta))

%{
    Now bound the derivative term.
%}


% controls

% Desired error bound to evaluate f and g.
desired_error_bound = 1e-16;


% bounds

pie = iv('pi');

% theta for stadium
theta = linspace(0,sup(2*pie),1001);
theta = iv(theta(1:end-1),theta(2:end));

% stadium for t_tilde
E_rho_t_tilde = (rho_t*exp(1i*theta)+exp(-1i*theta)/rho_t)/2;
% stadium for t
E_rho_t = (a_t+b_t)/2+(b_t-a_t)*E_rho_t_tilde/2;
% stadium for beta_tilde
E_rho_beta_tilde = (rho_beta*exp(1i*theta)+exp(-1i*theta)/rho_beta)/2;
% stadium for beta
E_rho_beta = (zf.a+zf.b)/2+(zf.b-zf.a)*E_rho_beta_tilde/2;

% get the maximum of 1-t on the stadium and check that is is less than 2.
max_E_rho_t_minus_one = max(abs(1-E_rho_t));
if sup(max_E_rho_t_minus_one) >= 2
   error('Theorem not valid'); 
end

% find the index k such that $|\beta| < k$ whenever $\beta$ is on or inside
% the stadium.
max_E_rho_beta = max(abs(E_rho_beta));

k_beta = ceil(sup(max_E_rho_beta));

%{
    Find a bound on f(t,beta) on the stadium
%}

M_f = iv(1); % first coefficient, $a_0 = 1$, of $f$.
an = iv(1);
for n = 1:k_beta
   an_bound = max(abs(an*((n^2+n-E_rho_beta)/(n+1)^2)/2));
   t_term_bound = max(abs((1-E_rho_t).^n));
   M_f = M_f + an_bound*t_term_bound; 
end
M_f = M_f + an_bound*max_E_rho_t_minus_one^(k_beta+1)/(2-max_E_rho_t_minus_one);

[N_t,N_beta,interpolation_error_t,interpolation_error_beta] = ...
    desired_N(rho_t,M_f,rho_beta,desired_error_bound);

interpolation_error_f = get_interpolation_error( ...
    N_beta,interpolation_error_beta,interpolation_error_t);

% Find a bound on f'(t,beta) on the stadiums
M_f_der = iv(1); % first coefficient, $a_0 = 1$, of $f$.
an = iv(1);
for n = 1:k_beta
   an_bound = max(abs(an*((n^2+n-E_rho_beta)/(n+1)^2)/2));
   t_term_bound = max(abs(n*(1-E_rho_t).^(n-1)));
   M_f_der = M_f_der + an_bound*t_term_bound*n; 
end

M_f_der = M_f_der + an_bound*(iv(1)/2)^(3-k_beta)/(2-max_E_rho_t_minus_one)^2;

[N_der_t,N_der_beta,interpolation_error_t,interpolation_error_beta] = ...
    desired_N(rho_t,M_f_der,rho_beta,desired_error_bound);

interpolation_error_f_der = get_interpolation_error( ...
    N_der_beta,interpolation_error_beta,interpolation_error_t);

x0_fun = @(beta) eval_standard(cf_standard,zf.a,zf.b,beta);

f0_fun = @(t,beta)f(t,beta,N_t,1,interpolation_error_f);
f1_fun = @(t,beta)f_der(t,beta,N_der_t,1,interpolation_error_f_der);

% sqrt_fun = @(beta) sqrt(1- x0_fun(beta).^2);

f0_der_fun = @(t,beta)f_der(t,beta,N_der_t,1,interpolation_error_f_der);
h0_fun_der= @(t,beta) f1_fun(t,beta)./(f0_der_fun(x0_fun(beta),beta).*sqrt(1-x0_fun(beta).^2));
h0_fun= @(t,beta) f0_fun(t,beta)./(f1_fun(x0_fun(beta),beta).*sqrt(1-x0_fun(beta).^2));
f0_grid = @(beta,t)f_grid(beta,t,N_t,iv(1),interpolation_error_f);
h0_grid = @(beta,t) h_grid(beta,t,f0_grid,f1_fun,x0_fun,interpolation_error_f);


