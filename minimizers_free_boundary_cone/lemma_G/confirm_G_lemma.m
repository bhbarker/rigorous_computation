function confirmed = confirm_G_lemma(G_der_poly,zeros_of_f_structure,a_beta,b_beta)

%{
    This part of the code evaluates a bound on the mininmum zero of $f$ on 
    the interval [a_beta,b_beta]. 
%}

% Break up the interval [a_beta,b_beta] into many small intervals.
dom = linspace(a_beta,b_beta,1001);
dom = iv(dom(1:end-1),dom(2:end));

% We evaluate the interpolation of the zeros of f rigorously.
dom_tilde = (2/(zeros_of_f_structure.b-zeros_of_f_structure.a))* ...
    (dom-(zeros_of_f_structure.a+zeros_of_f_structure.b)/2);
theta_dom = acos(dom_tilde);
ind = 0:1:zeros_of_f_structure.N-1;
T = cos(theta_dom.'*ind);
f = T*zeros_of_f_structure.cf_p;

% We get a lower bound on the zeros of f on the interval [a_beta,b_beta]
a_t = min(inf(f-zeros_of_f_structure.root_error_bound));

% We verify that g^3(t,beta)G'(t,beta) > 0 for t in [a_t,1] and beta in [a_beta,b_beta].

b_t = 1;
nodes_beta = 100;
nodes_t = 1000;
confirmed = G_der_is_positive(G_der_poly,a_beta,b_beta,nodes_beta,a_t,b_t,nodes_t);










