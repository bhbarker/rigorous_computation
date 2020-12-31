function out = eval_standard(cf,a,b,x)

x = (2*x-(a+b))/(b-a);


out = cf(end-1)+cf(end)*x;

for j = 2:length(cf)-1
    out = cf(end-j)+out.*x;
end

