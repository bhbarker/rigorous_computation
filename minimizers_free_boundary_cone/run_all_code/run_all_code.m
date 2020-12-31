% fun all code for the proofs
clc; close all; clear all; curr_dir = cd;

%{
    You must first start intlab.
%}


% set the paths for running the code
cd('..');
set_paths;

% Find the zeros of f
cd('bound_zero_of_f');
% bound_zero_of_f

% verify Lemma G
cd('../lemma_G');
% lemma_G;

% very super solutions
cd('../super_solution');
driver_batch_jobs




cd(curr_dir);