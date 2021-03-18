function [cf,fun] = get_chebyshev_coefficients(f,a,b,opt)
% =================================================================
%Purpose of the function:
%get the chebyshev_coefficients of function f
% Parameters:
% f: evaluated points of f
% a: left end point of the interval
% b: right end point of the interval
% opt: option of methods of interpolation
% =================================================================
sz = size(f);
if length(sz)>2
   error('input must be a vector'); 
end
if sz(2)>sz(1)
    f = f.';
end
if size(f,2) > 1
   error('input must be a vector');
end

if strcmp(opt,'first kind')
    
    N = length(f);
    Id2 = (2/N)*speye(N);
    Id2(1,1) = Id2(1,1)/2;
    theta = ((0:1:N-1)+0.5)*pi/N;
    Tcf = Id2*cos(theta.'*(0:1:N-1)).';
    cf = Tcf*f;
    T_x = (2/(b-a))*(repmat((0:1:N-1),N,1) ...
    .*sin(theta.'*(0:1:N-1)))./sin(theta.'*ones(1,N));
    fx = T_x*cf;
    cf_x = Tcf*fx;
    fun = @(x)eval_cf(x,cf,cf_x,a,b);
    
elseif strcmp(opt,'FirstKindExtrema')
    
    
    
end


% Transformation to get Chebyshev coefficients

function [out1,out2] = eval_cf(x,cf,cf_x,a_x,b_x)


    N = length(cf);
    xtilde = (x-0.5*(a_x+b_x))/(0.5*(b_x-a_x));
    theta = acos(xtilde);
    T = cos(theta.'*(0:1:N-1));
    out1 = T*cf;
    out2 = (T*cf_x).';
    
    
   







