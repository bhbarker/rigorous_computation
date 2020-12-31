clc; close all; beep off; clear all; commandwindow; curr_dir = cd;

intvalinit('displayinfsup');
% intvalinit('displaymidrad');


d{1}.left = '0.41';
d{1}.right = '0.42';


% coefficients for M
zk_data.zk_standard = [iv('0.193');iv('1.005');iv('-0.09');iv('0.03')];
zk_data.mu = iv('0.1');
zk_data.x0 = [ -0.214618461719847; 0.653405380267736];


% Load the coefficients $z_k$.
curr_dir = cd;
cd('../data');
% ld = load(file_name);
% zk_data = ld.data; % coefficients $z_k$
ld = load('data_bound_zero_of_f'); % load zeros of f data
zf = ld.data;
cd(curr_dir);


for j = 1:length(d)
    j
        
        % interval of c values
        c_int = [iv(d{j}.left),iv(d{j}.right)];

        % file name containing the coefficients
        file_name = ['c_from_',d{j}.left,'_to_',d{j}.right]; 
        file_name = strrep(file_name,'.','p');


        small_interval = iv(-1e-3,1e-3);

        % -----------------------

        verify_super_solution(file_name,c_int,small_interval,zk_data,zf);
    
end





