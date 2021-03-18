function out = iv(num1,num2)
% =================================================================
%Purpose of the function:
% returns a mathematically rigorous enclosure of the number of interval given.
% Parameters:
% num1: lower bound
% num2: upper bound
% =================================================================


if nargin == 1
    out = intval(num1);
elseif nargin == 2
    out = hull(intval(num1),intval(num2));
end