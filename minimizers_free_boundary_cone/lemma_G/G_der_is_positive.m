function success = G_der_is_positive(G_der_poly,a_beta,b_beta,nodes_beta,a_t,b_t,nodes_t)
% =================================================================
%Purpose of the function:
% verifies the inequality given by (4.6)

% Parameters:
% G_der_poly: the class that contains information about rigorous
% calculation of G'
% a_beta: lower bound of beta
% b_beta: upper bound of beta
% nodes_beta: number of nodes for interpolation with respect to beta
% a_beta: lower bound of t
% b_beta: upper bound of t
% nodes_t: number of nodes for interpolation with respect to t
% =================================================================

% create a subgrid of intrevals on which to verify sign of g^3G'
y = linspace(a_beta,b_beta,nodes_beta);
x = linspace(a_t,b_t,nodes_t);

x = iv(x(1:end-1),x(2:end));
y = iv(y(1:end-1),y(2:end));

% Evaluate interpolating polynomial of g^3G' rigorously on the grid found.
M = G_der_poly.grid_eval(x,y);

% Get the minimum of the rigorous interpolation
Z = (inf(M));

% Determine if the rigorous interpolation is less than zero.
ind = find(Z<0);

% check if g^3G' was always positive
if isempty(ind)
    % set confirmation of verifying the result to true
   success = 1; 
else
    % set confirmation of verifying the result to false
    success = -1;
end



