function R = bound_R_for_M_zero(beta,zk)

% make sure that zk is a column vector
if size(zk,2)>1
   zk = zk.'; 
end

% make sure that beta is a column vector
if size(beta,2) > 1
   beta = beta.'; 
end


N = length(zk); % the number of coefficients, the number of Legndre polys
k = iv(1:1:N-1); % index (not including zero index)

bk = (-1 + sqrt(1+(1./beta)*(8*(k+1).*k)))/2; % The powers of R

R = iv(1);
other_terms_bound = iv(max(sup(R.^bk)))*zk(2:end);
while sup(other_terms_bound) >= inf(abs(zk(1)))
    R = iv(sup(R/2));
    other_terms_bound = iv(max(sup(R.^bk)))*zk(2:end);
end








