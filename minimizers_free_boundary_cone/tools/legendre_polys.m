function cf = legendre_polys(N,make_rigorous)
% =================================================================
%Purpose of the function:
%   each row corresponds to the coefficints of the Legendre polynomial 
%   in the standard basis centered at x = 0.
%     
%   For example, legendre_polys(3) returns [[1,0,0],[0,1,0],[-0.5,0, 1.5]]

% Parameters:
% N: number of legendre polynomial
% make_rigorous: option for double or interval arithmetic
% =================================================================

    if strcmp(make_rigorous,'on')
        cf = iv(zeros(N,N));
        cf(1,1) = iv(1); % first Legendre polynomial
        cf(2,2) = iv(1); % second Legendre polynomial
        two = iv(2);
        zero = iv(0);
    else
        cf = zeros(N,N);
        cf(1,1) = 1;
        cf(2,2) = 1;
        two = 2;
        zero = 0;
    end

    
    % obtain the other Legendre polynomials
    for n = 1:N-2
        
        ind = n+1;
       
        % shift the coefficints of the nth Chebyshev polynomial one to the right
        r = [zero,cf(ind,1:end-1)];
        
        % The (n+1)th Legendre polynomial
        cf(ind+1,:) = ((two*n+1)/(n+1))*r-(n*cf(ind-1,:))/(n+1);
        
    end







