function interpolation_error = get_interpolation_error(N,error_beta,error_t)

pie = iv('pi');
LAM = (2/pie)*(log(N)+iv('0.58')+log(8/pie))+pie/(72*N^2);
interpolation_error = error_beta+LAM*error_t;