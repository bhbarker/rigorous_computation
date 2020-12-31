commandwindow


% Prove that G_c'(t) > 0 for t_c <= t < 1
clc; close all; beep off; clear all
% intvalinit('DisplayMidRad');
intvalinit('displayinfsup');
clc;

%% controls

% Interval of which t is an element.

a_t = iv('-0.204'); % left end point of t interval
b_t = iv(1); % right end point of t interval

a_beta = iv('-0.161');% left end point of beta interval
b_beta = iv(2);  % right end point of beta interval

% Desired error bound to evaluate f and g.
desired_error_bound = 1e-16;

% rho giving stadium
rho_t = iv('2.9');
rho_beta = iv(30);


%% bounds
pie = iv('pi');

% t interval
t_interval = iv(a_t,b_t);

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

max_E_rho_t = max(abs(E_rho_t));

% find the index k such that $|\beta| < k$ whenever $\beta$ is on or inside
% the stadium.
max_E_rho_beta = max(abs(E_rho_beta));


k_beta = ceil(sup(max_E_rho_beta));

%{
    Find a bound on f(t,beta) on the stadiums
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

%{
    Now bound the derivative term.
%}

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

%{
    Now bound the second derivative term.
%}

% Find a bound on f''(t,beta) on the stadiums
M_f_der_der = iv(1); % first coefficient, $a_0 = 1$, of $f$.
an = iv(1);
for n = 1:k_beta
   an_bound = max(abs(an*((n^2+n-E_rho_beta)/(n+1)^2)/2));
   t_term_bound = max(abs(n*(n-1)*(1-E_rho_t).^(n-2)));
   M_f_der_der = M_f_der_der + an_bound*t_term_bound*n*(n-1); 
end

M_f_der_der = M_f_der_der + an_bound*(iv(1)/2)^(6-k_beta)/(2-max_E_rho_t_minus_one)^3;


[N_der_der_t,N_der_der_beta,interpolation_error_t,interpolation_error_beta] = ...
    desired_N(rho_t,M_f_der,rho_beta,desired_error_bound);

interpolation_error_f_der_der = get_interpolation_error( ...
    N_der_der_beta,interpolation_error_beta,interpolation_error_t);


%{
    Get a bound on G(t,beta)
%}

M1 = M_f;
M2 = M_f_der;
M3 = M_f_der_der;

sigma_bound = (iv(9)/8)*max_E_rho_beta;
prod1 = 2*sigma_bound*M1*M2*M1^3;
prod2 = 8*max_E_rho_t*(M1*M2)^2*M1;
prod3 = 4*(1+max_E_rho_t^2)*(M1*M2)*(M3*M1^2+(M1*M2)*M1+M1^2*M3+M1*M2^2);

M_G = prod1+prod2+prod3;


N_G_t = 1;
eta_t = log(rho_t);
D_rho_t = (rho_t+1/rho_t)/2-1;
L_rho_t = sqrt(rho_t^2+1/rho_t^2);
constant = M_G*L_rho_t/(D_rho_t);
interpolation_error_t = iv(N_G_t)+1;
while sup(interpolation_error_t) >= desired_error_bound
    N_G_t = N_G_t + 1;
    interpolation_error_t = constant/sinh(eta_t*(N_G_t+1));
end

N_G_beta = 1;
eta_beta = log(rho_beta);
D_rho_beta = (rho_beta+1/rho_beta)/2-1;
L_rho_beta = sqrt(rho_beta^2+1/rho_beta^2);
constant = M_G*L_rho_beta/(D_rho_beta);
interpolation_error_beta = iv(N_G_beta)+1;
while sup(interpolation_error_beta) >= desired_error_bound
    N_G_beta = N_G_beta + 1;  
    interpolation_error_beta = constant/sinh(eta_beta*(N_G_beta+1));
end

LAM = (2/pie)*(log(N_G_beta)+iv('0.58')+log(8/pie))+pie/(72*N_G_beta^2);
interpolation_error_G = interpolation_error_beta+LAM*interpolation_error_t;

%{
    Form G(t,beta) and interpolate it.
%}

f0_fun = @(t,beta)f_fun(t,beta,N_t,iv(1),interpolation_error_f);
f1_fun = @(t,beta)f_der_fun(t,beta,N_der_t,iv(1),interpolation_error_f_der);
f2_fun = @(t,beta)f_der_der_fun(t,beta,N_der_der_t,iv(1),interpolation_error_f_der_der);

g0_fun = @(t,beta)f_fun(t,-beta/8,N_t,iv(1),interpolation_error_f);
g1_fun = @(t,beta)f_der_fun(t,-beta/8,N_der_t,iv(1),interpolation_error_f_der);
g2_fun = @(t,beta)f_der_der_fun(t,-beta/8,N_der_der_t,iv(1),interpolation_error_f_der_der);


G_der_fun = @(t,beta)G_der_g3(t,beta,f0_fun,f1_fun,f2_fun,g0_fun,g1_fun,g2_fun);



t = linspace(0.69,0.7,11);
t = iv(t(1:end-1),t(2:end));
beta = linspace(1.6,1.61,11);
beta = iv(beta(1:end-1),beta(2:end));
temp = G_der_fun(t,beta).';
% left = inf(temp(1,1));
% right = sup(temp(1,1));

G_der_poly = iv_poly([a_t;a_beta],[b_t;b_beta],G_der_fun,interpolation_error_G, ...
    N_G_t,N_G_beta);

curr_dir = cd;
cd('../data');
ld = load('data_bound_zero_of_f');
zeros_of_f_structure = ld.data;
zeros_of_f_structure.N = length(zeros_of_f_structure.cf_p);
cd(curr_dir);

R_beta_point = '0.58828';
beta_vals = [inf(2/(1+iv('0.5882')^2)), sup(2/(1+iv(R_beta_point)^2));
inf(2/(1+iv('0.588')^2)), sup(2/(1+iv('0.5882')^2));
inf(2/(1+iv('0.58')^2)), sup(2/(1+iv('0.588')^2));
inf(2/(1+iv('0.5')^2)), sup(2/(1+iv('0.58')^2));
inf(2/(1+iv('0')^2)), sup(2/(1+iv('0.5')^2)) ];

for j = 1:size(beta_vals,1)
    
    confirmed = confirm_G_lemma(G_der_poly,zeros_of_f_structure,...
        beta_vals(j,1),beta_vals(j,2));

    if ~(confirmed == 1)
       error('The lemma is not established.'); 
    end
end

fprintf(['Lemma G has been verified for beta in [0,',R_beta_point,'].\n\n']);







