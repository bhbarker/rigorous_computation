function out = eval_standard(cf,a,b,x)
% =================================================================
%Purpose of the function:
%evaluate function based on interpolation coefficients in standard basis
% Parameters:
% cf: the standard basis coefficients
% a: left end point of the interval
% b: right end point of the interval
% x: values to evaluate
% =================================================================

%shift the interval to {-1,1]

x = (2*x-(a+b))/(b-a);


out = cf(end-1)+cf(end)*x;

for j = 2:length(cf)-1
    out = cf(end-j)+out.*x;
end

