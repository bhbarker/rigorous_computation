function [N_t,N_beta,interpolation_error_t,interpolation_error_beta] = desired_N(rho_t,M,rho_beta,desired_error_bound)


N_t = 1;
eta_t = log(rho_t);
D_rho_t = (rho_t+1/rho_t)/2-1;
L_rho_t = sqrt(rho_t^2+1/rho_t^2);
constant = M*L_rho_t/(D_rho_t);
interpolation_error_beta = iv(N_t);
while sup(interpolation_error_beta) >= desired_error_bound
    N_t = N_t + 1;
    interpolation_error_beta = constant/sinh(eta_t*(N_t+1));
end


N_beta = 1;
eta_beta = log(rho_beta);
D_rho_beta = (rho_beta+1/rho_beta)/2-1;
L_rho_beta = sqrt(rho_beta^2+1/rho_beta^2);
constant = M*L_rho_beta/(D_rho_beta);
interpolation_error_t = iv(N_beta);
while sup(interpolation_error_t) >= desired_error_bound
    N_beta = N_beta + 1;  
    interpolation_error_t = constant/sinh(eta_beta*(N_beta+1));
end




