close all; clc; clear all; curr_dir = cd;
intvalinit('displayinfsup')


%----------------------------------------------------------
% Approximate the curve of zeros of f, that is f(x(c)) = 0
%----------------------------------------------------------

% The number of terms in series representation of f using double arithmetic
num_terms = 80;

a = 1.47; % left value of beta
b = 2; % right value of beta
N = 13; % number of Chebyshev nodes to inerpolate zero of f 

plot_it = 'on';

%
% use double arithmetic to get a polynomial, p, approximating zeros of f
%

theta = ((0:1:N-1)+0.5)*pi/N;
dom = (a+b)/2+(b-a)*cos(theta)/2; % interpolation nodes in original interval 

zr = zeros(length(theta),1); % zeros of f at beta nodes
iter = 40; % number of iterations in the Newton solver to find zeros of f
for j = 1:length(theta)
    beta = dom(j); 
    fun = @(t)f_newton(t,beta,num_terms,1);
    zr(j) = newton(fun,0,iter); % get zero of f.
end
[cf_double,fun] = get_chebyshev_coefficients(zr,a,b,'first kind');

if strcmp(plot_it,'on')
    % plot the approximation of the zeros of f
    figure;
    hold on;
    plot(dom,zr,'-k','LineWidth',2);
    plot(dom,fun(dom),'.r','MarkerSize',18);
    h = xlabel('\beta');
    set(h,'FontSize',18);
    h = ylabel('x_0(\beta)');
    set(h,'FontSize',18);
    h = gca;
    set(h,'FontSize',18);
end


%-------------------------------------------------------- 
% Make bound rigorous
%-------------------------------------------------------- 

% turn coefficients into intervals
cf_p = iv(cf_double);
a = iv(a);
b = iv(b);

% interval constants
one = iv(1);
two = iv(2);
half = one/2;
pie = iv('pi');



%{
    Assumption 3 of Lemma. 
    Verify that $a\leq p(\beta) -\delta$ and $p(\beta)+\delta \leq b$.
%}

% Get upper bounds on p(beta), p(beta) approximates the zero of f
x = linspace(a.inf,b.sup,1000); % get mesh on interval [a,b]
x = iv(x(1:end-1),x(2:end)); % make into intervals
xtilde = (2*x-a-b)/(b-a); % coordinate change for polynomial evaluation
theta_tilde = acos(xtilde); % coordinate change 
p_tilde = cos(theta_tilde.'*(0:1:N-1))*cf_p; % evaluate p(beta)
p_upper = max(p_tilde); % find upper bound on p(beta)
p_lower = min(p_tilde); % find lower bound on p(beta)

% A bound on how large delta can be.
delta_max = iv(1e-6);

% find upper and lower bound on t when $|t-p(beta)|\leq delta_max$
t_upper = p_upper+delta_max;
t_lower = p_lower-delta_max;

% check the assumption of the lemma that t_lower >= a, t_upper <= b
% That is, verify that $a\leq p(\beta) -\delta$ and $p(\beta)+\delta \leq b$,
if (t_upper.sup > 3) == 1
    error('t_upper <= 1 required');
end

if (t_lower.inf < -1) == 1
    error('t_lower >= -1 required');
end


%{
    Assumption 1 of the Lemma. Verify that
    $|f(p(\beta),\beta)|<\varepsilon$.
%}

% bound on the t values we consider, i.e. $|1-t|< r < 2$.
r = max(abs(1-t_upper),abs(1-t_lower));

% remainder error for series representation of f
f_err = 2*(r/2)^num_terms/(1-(r/2));
f_err = iv(-f_err,f_err);

% degree of polynomial
N_x = N*num_terms;

% new points
theta_x = ((0:1:N_x-1)+half)*pie/N_x;
dom_x = (a+b)/2+(b-a)*cos(theta_x)/2;

zp = cos(theta_x.'*(0:1:N-1))*cf_p;

% get values of f
fv = f(zp.',dom_x,num_terms,one);
fv = fv+f_err;
eps = iv(sup(max(abs(fv)))); % Verify that $|f(p(\beta),\beta)|<\varepsilon$

%{
    Assumption 2 of the lemma. 

    Show that $|\frac{\partial f}{\partial t}(s,\beta)|>k$ for all $s\in\R$ 
    such that $|s-p(\beta)|<\delta$ where $\delta:= \varepsilon/k$,

%} 

if t_upper.sup >= 0
   error('The term below is not correct.'); 
end

%  This bounds the derivative of $f$ on the domain of interest since 
% $a_n(\beta) <= 0$ for $n\geq 1$, and since $t<= 0$. 
k_der_bound = a/2; 

delta = eps/k_der_bound;

if delta.sup > delta_max.inf
    error('delta > delta_max');
end

% save the results
data.cf_p = cf_p;
data.root_error_bound = delta;
data.a = a;
data.b = b;
cd('../data');
save('data_bound_zero_of_f','data');
cd(curr_dir);


fprintf('\n\nBound zeros of f Lemma was successful.\n\n')



