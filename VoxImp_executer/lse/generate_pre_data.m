clc; close all; clear all; format long e;

% -------------------------------------------------------------------------
%                  Add the Current Path to Workspace
% -------------------------------------------------------------------------
parpool('local',24);
pre_define_the_path_for_folders

dx=1;
ko=1;
fl_no_fft=1;



L=9000; M=9000; N=90;
[T_cap]=prestore_T_cap(dx,L,M,N);
[T_henry] = prestore_T_henry(ko,L,M,N,1e-8);

save('prestore_9000_9000_90.mat','T_cap','T_henry')


% L=3000; M=2000; N=10;
% [T_cap]=prestore_T_cap(dx,L,M,N);
% [T_henry] = prestore_T_henry(ko,L,M,N,1e-8);
% 
% save('prestore_3000_2000_10.mat','T_cap','T_henry')
