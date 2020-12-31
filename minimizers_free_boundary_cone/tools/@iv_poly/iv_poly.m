classdef iv_poly
    
    properties
        a % left end points
        b % right end points
        cf % coefficients
        N1 % degree+1 in first variable
        N2 % degree+1 in second variable
    end
    
    methods
        
        % instantiation
        function obj = iv_poly(a,b,cf_or_fun,interpolation_error,N1,N2)
            
            if nargin > 3
                % get the Chebyshev coefficients
                cf = get_coefficients(a,b,cf_or_fun,N1,N2);
                cf(1,1) = cf(1,1)+interpolation_error; 
            else
                % Polynomial defined with Chebyshev coefficients
                cf = cf_or_fun;
            end
            
            obj.a = a;
            obj.b = b;
            obj.cf = cf;
            obj.N1 = size(cf,1);
            obj.N2 = size(cf,2);
        end
        
        % addition of two polynomials
        function out = plus(p1,p2)
    
            if ~(p1.a == p2.a)
                error('Polynomials do not share left end points');
            elseif ~(p1.b == p2.b)
                error('Polynomials do not share right end points');
            end
    
            out = iv_poly(p1.a,p1.b,cheby_add(p1.cf,p2.cf));
            
        end

        % addition of two polynomials
        function out = mtimes(p1,p2)
    
            if ~(p1.a == p2.a)
                error('Polynomials do not share left end points');
            elseif ~(p1.b == p2.b)
                error('Polynomials do not share right end points');
            end
   
            out = iv_poly(p1.a,p1.b,cheby_product(p1,p2));
            
        end
        
        % evaluation of a polynomial
        function out = grid_eval(p,x,y)
           
            if size(x,1) > 1
                x = x.';
            end
            if size(y,1) > 1
               y = y.'; 
            end

%             if x(1) < p.a(1)
%                 error('Evaluating outside the valid domain.');
%             end
%             if x(end) > p.b(1)
%                 error('Evaluating outside the valid domain.');
%             end
%             if y(1) < p.a(2)
%                 error('Evaluating outside the valid domain.');
%             end
%             if y(end) > p.b(2)
%                 error('Evaluating outside the valid domain.');
%             end
            
            xtilde = intersect((2*x-(p.a(1)+p.b(1)))/(p.b(1)-p.a(1)),iv(-1,1));
            theta_x = acos(xtilde);

            ytilde = intersect((2*y-(p.a(2)+p.b(2)))/(p.b(2)-p.a(2)),iv(-1,1));
            theta_y = acos(ytilde);

            % 2d Chebyshev polynomials
            T_x = (cos((0:1:p.N1-1).'*theta_x)).'; 
            T_y = (cos(theta_y.'*(0:1:p.N2-1))).';
            
            out = T_x*p.cf*T_y;
            
        end
        

        
    end
end




%             obj.a = a;
%             obj.b = b;
%             obj.cf = cf;
%             obj.N1 = size(cf,1);
%             obj.N2 = size(cf,2);
% ---------------------------------------------------
% get_coefficients
% ---------------------------------------------------

function cf = get_coefficients(a,b,fun,N1,N2)


    half = 1/iv(2);
    pie = iv('pi');
    theta1 = ((0:1:(N1-1))+half)*pie/N1;
    theta2 = ((0:1:(N2-1))+half)*pie/N2;

    % get domain points
    xtilde = cos(theta1);
    x = half*(a(1)+b(1))+half*(b(1)-a(1))*xtilde;
    ytilde = cos(theta2);
    y = half*(a(2)+b(2))+half*(b(2)-a(2))*ytilde;
    
    % 2d Chebyshev polynomials
    T1_cf = cos((0:1:N1-1).'*theta1); 
    T2_cf = cos(theta2.'*(0:1:N2-1));
    
    
    cf= (4/(N1*N2))*T1_cf*fun(x,y)*T2_cf;
    cf(:,1) = cf(:,1)/2;
    cf(1,:) = cf(1,:)/2;
    
end

% ---------------------------------------------------
% Multiply two Chebyshev polynomials together
% ---------------------------------------------------

function cf = cheby_product(f,g)

    % degere of the new polynomial
    deg_x = f.N1+g.N1-2;
    N1 = deg_x - 1;
    deg_y = f.N2+g.N2-2;
    N2 = deg_y - 1;
    
    fun = @(x,y) f.grid_eval(x,y).*g.grid_eval(x,y);
    cf = get_coefficients(f.a,f.b,fun,N1,N2);

end

% ---------------------------------------------------
% Add two Chebyshev polynomials together
% ---------------------------------------------------

function cf = cheby_add(f,g)
% add two Chebyshev polynomials with the same domain

    [sfx,sfy] = size(f);
    [sgx,sgy] = size(g);

    m = max(sfx,sgx);
    n = max(sfy,sgy);

    if m > sfx
        f = [f; iv(zeros(m-sfx,sfy))];
    end
    if n > sfy
        f = [f, iv(zeros(m,n-sfy))];
    end

    if m > sgx
        g = [g; iv(zeros(m-sgx,sgy))];
    end
    if n > sgy
        g = [g, iv(zeros(m,n-sgy))];
    end

    cf = f+g;

end

% ---------------------------------------------------
% Local copy needed
% ---------------------------------------------------
function out = iv(num1,num2)
% function out = iv(num)
%
% returns a mathematically rigorous enclosure of the number of interval
% given.

    if nargin == 1
        out = intval(num1);
    elseif nargin == 2
        out = hull(intval(num1),intval(num2));
    end

end




