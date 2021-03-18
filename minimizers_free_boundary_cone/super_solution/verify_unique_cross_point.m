function verify_unique_cross_point(R,Y,beta,zk_standard,cf_legendre_standard_basis)
% =================================================================
%Purpose of the function:
%verify the inequality given in Lemma 5.1 (3)
% Parameters:
% file_name (string): Data is saved as "verified_" concatenated with file_name.
% c_int (2x1 array of intervals): Contains two intervals, the left and right
% end point of the interval of c-values that will be verified.
% small_interval (interval): The cross point will be rigorously verified by
% initializing interval Newton's method with an approximation of the cross
% point plus the small_interval. 
% zk_data (array of intervals): The coefficients for M.
% zf (structure): The structure that contains data about the zeros of f.
% =================================================================

iv_fun = @(x)iv(x);

% return
%{
 Y is a vector containing the y-coordinate of the cross point corresponding
 to each value of $\beta$ in the vector beta. We call y0 the max of these
 y-values.
%}
y0 = max(sup(Y));

%{
    Now we split up the interval [0,y0] into evenly spaced subintervals. We
    will show that for each subinterval, if M(x,y) = 0, then the gradient
    condition is satisfied.
%}
num_intervals_y = 30;
yv = linspace(0,y0,num_intervals_y+1);
yvals = iv(yv(1:end-1),yv(2:end));

tic

for j = 1:length(yvals)
    
    

    % make y a constant vector the size of beta for input into M_fun
    y = yvals(j)*ones(size(beta));
    
    % number of intervals in the x variable
    num_intervals_x = 30;

    % get the admissible values of x
    xv = linspace(-1,-sup(sqrt(iv(max(0,inf(R^2-y(j)^2)),num_intervals_x+1))));
    xvals = iv(xv(1:end-1),xv(2:end));

    %{
        Iterate through each (x,y) interval pairs to verify condition. Note
        that beta is a large vector, so vectorization is occurring.
    %}
    for k = 1:length(xvals)


        % make x a constant vector the size of beta for input into M_fun
        x = xvals(k)*ones(size(beta));

        % compute M(x,y,beta)
        M = M_fun(x,y,beta,zk_standard,cf_legendre_standard_basis,iv_fun);

        if max(isinf(M)) == 1
            error('Failed to verify');
        end

        if max(isnan(M)) == 1
           error('Failed to verify'); 
        end

        % find where M might have a zero
        ind = find(sup(M).*inf(M) <= 0);

        % Where M might have a zero, check the gradient condition is
        % satisfied.

        if ~isempty(ind)

            grad_M_val = grad_M_fun(x(ind),y(ind),beta(ind),zk_standard,cf_legendre_standard_basis,iv_fun);

            if max(isnan(grad_M_val)) == 1
               error('Failed to verify'); 
            end
            if max(isinf(grad_M_val)) == 1
                error('Failed to verify');
            end
            if max(sup(abs(grad_M_val))>= 1) == 1
                error('Failed to verify');
            end

        end

    end
            
end

toc