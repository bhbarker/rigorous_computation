function [g0_fun,g0_fun_der,g0_grid] = interp_g(a_beta,b_beta,a_t,b_t,rho_beta,rho_t)

% =================================================================
%Purpose of the function:
%interpolation function g and g' with desired error
% Parameters:
% a_beta: lower bound of beta
% b_beta: upper bound of beta
% a_t: lower bound of t
% b_t: upper bound of t
% rho_t: interpolation variable
% rho_beta: interpolation variable
% =================================================================

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
E_rho_beta = (a_beta+b_beta)/2+(b_beta-a_beta)*E_rho_beta_tilde/2;

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




f1_fun = @(t,beta)f_der(t,beta,N_der_t,1,interpolation_error_f_der);
g0_fun = @(t,beta)f(t,-beta/8,N_t,1,interpolation_error_f);
g0_fun_der = @(t,beta)f1_fun(t,-beta/8);
f0_grid = @(beta,t)f_grid(beta,t,N_t,iv(1),interpolation_error_f);
g0_grid = @(beta,t) f0_grid(-beta/8,t);


