% use recent DIRECTFN (faster)
use_recent_DIRECTFN=0;

% find the current folder
currentfolder = pwd;

% get system-dependent path separator
psep = pathsep();

% include all the needed path (addpath() adds recursively),
% but do NOT include the Git directory


p1  = genpath([currentfolder, filesep, 'SVD_Tuck']);
p2  = genpath([currentfolder, filesep, 'src_lin_vie']);
p3  = genpath([currentfolder, filesep, 'src_iterative_solver']);
p4  = genpath([currentfolder, filesep, 'pre_process']);
p5  = genpath([currentfolder, filesep, 'one_over_R_Cap']);
p6  = genpath([currentfolder, filesep, 'lse']);
p8  = genpath([currentfolder, filesep, 'STRUMPACK']);
addpath([p1, psep,p2, psep, p3, psep, p4, psep, p5, psep, p6, psep, p8, psep ]);

if (ispc)
  p7 = genpath([currentfolder, filesep, 'windows_mex']);
else
  p7 = genpath([currentfolder, filesep, 'Linux_mex']);
end
addpath(p7);

% get system-dependent file separator
filesep = filesep();

currentfolder = pwd;