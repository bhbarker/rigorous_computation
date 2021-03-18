function out = standard_basis_from_chebyshev(cf,make_rigorous)
% =================================================================
%Purpose of the function:
% find the coefficients on the standard basis given coefficients on the
% chebshev basis
% Parameters: 
% cf: coefficients on the chebshev basis
% make_rigorous: option to use interval arithmetic
% =================================================================
N = length(cf);
if size(cf,2) > size(cf,1)
   cf = cf.'; 
end

if strcmp(make_rigorous,'on')
    out = iv(zeros(N));
    out(1,1) = iv(1);
    out(2,2) = iv(1);
    two = iv(2);
else
    out = zeros(N);
    out(1,1) = 1;
    out(2,2) = 1; 
    two = 2;
end

for j = 3:N
   out(j,:) = two*[0,out(j-1,1:N-1)]-out(j-2,:); %T_{n+1}= 2*x*T_n-T_{n-1}
end
out = sum(diag(cf)*out);
