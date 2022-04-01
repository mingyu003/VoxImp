function [idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams_Imp(r,dx,Cnt,Dims,Orients,e_r,s_e,fl_plot_intcons)
%%    Defines Constitutive parameters of interconnects
% _________________________________________________________________________
%
%       Defines interconnects in a given domain.
%       Optionally assigns constitutive parameters to interconnects
%
% _________________________________________________________________________
%
%% INPUT
%   r           4D (LxMxNx3) array with domain voxelized grid coordinates
%   Cnt         Cartesian coordinates of the centers of interconnects
%               row: interconnect id, column: dimensions
%   Dims        Dimensions of interconnects (LxWxH) (along x, y, and z)
%               row: interconnect id, column: length,width,height
%   Orients     Orientations of conductors, string
%   e_r         value of relative epsilon
%   s_e         value of electric conductivity
%
%% OUTPUT
%   idx         indices of the positions, 2nd dimension indicates which
%               conductor the voxel is located in
%   epsilon_r   3D (LxMxN) array with relative epsilon
%   sigma_e     3D (LxMxN) array with electric conductivity
%   grid_intcon 4D (LxMxNx5), 4-dim the first 3 entries mean coordinates,
%   the 4th entry means the numbering of conductor, 5th entry is the diel
%   numbering
% -------------------------------------------------------------------------
%
%   Mingyu Wang -- mingyu003@e.ntu.edu.sg
%   Computational Electromagnetics Group, NTU
%
% _________________________________________________________________________


% -------------------------------------------------------------------------
% Initialization of variables
% -------------------------------------------------------------------------

if(nargin < 4 )
    fprintf(1, '\n ERROR: not enough arguments\n');
    return
end
if(nargin < 5 || isempty(e_r))
    e_r = 1;
end
if(nargin < 6 || isempty(s_e))
    s_e = 0;
end

% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------
num_conds=max(Cnt(:,4));
[L,M,N,~] = size(r);

num_intcon=size(Cnt,1);

% define bounds for each interconnect
x_bnd=zeros(num_intcon,4);y_bnd=zeros(num_intcon,4);z_bnd=zeros(num_intcon,4);%3rd col: cond numbering, 4th col: diel numbering
for kk=1:num_intcon
    temp_cen=Cnt(kk,1:3);
    temp_dim=Dims(kk,1:3);
    cond_ind=Cnt(kk,4);
    diel_ind=Cnt(kk,5);
    if (Orients(kk,1) == 'x')
        x_bnd(kk,1:4)=[temp_cen(1)-temp_dim(1)*0.5 temp_cen(1)+temp_dim(1)*0.5 cond_ind diel_ind];
        y_bnd(kk,1:4)=[temp_cen(2)-temp_dim(2)*0.5 temp_cen(2)+temp_dim(2)*0.5 cond_ind diel_ind];
        z_bnd(kk,1:4)=[temp_cen(3)-temp_dim(3)*0.5 temp_cen(3)+temp_dim(3)*0.5 cond_ind diel_ind];
    elseif (Orients(kk,1) == 'y')
        x_bnd(kk,1:4)=[temp_cen(1)-temp_dim(2)*0.5 temp_cen(1)+temp_dim(2)*0.5 cond_ind diel_ind]; 
        y_bnd(kk,1:4)=[temp_cen(2)-temp_dim(1)*0.5 temp_cen(2)+temp_dim(1)*0.5 cond_ind diel_ind];
        z_bnd(kk,1:4)=[temp_cen(3)-temp_dim(3)*0.5 temp_cen(3)+temp_dim(3)*0.5 cond_ind diel_ind];
    elseif (Orients(kk,1) == 'z')
        x_bnd(kk,1:4)=[temp_cen(1)-temp_dim(3)*0.5 temp_cen(1)+temp_dim(3)*0.5 cond_ind diel_ind];
        y_bnd(kk,1:4)=[temp_cen(2)-temp_dim(2)*0.5 temp_cen(2)+temp_dim(2)*0.5 cond_ind diel_ind];
        z_bnd(kk,1:4)=[temp_cen(3)-temp_dim(1)*0.5 temp_cen(3)+temp_dim(1)*0.5 cond_ind diel_ind];
    else
        error('Orients should have only x, y, or z strings')
    end
    
end


tic;
tola=1e-12;
idx_temp=cell(num_intcon,1);
for kk=1:num_intcon
    boolean_tens=zeros(L*M*N,1); 
    for ll=round((x_bnd(kk,1))/dx+1):round((x_bnd(kk,2))/dx)
        for mm=round((y_bnd(kk,1))/dx+1):round((y_bnd(kk,2))/dx)
            for nn=round((z_bnd(kk,1))/dx+1):round((z_bnd(kk,2))/dx)
                if ( r(ll,mm,nn,1) > x_bnd(kk,1)-tola && r(ll,mm,nn,1) < x_bnd(kk,2)+tola &&...
                        r(ll,mm,nn,2) > y_bnd(kk,1)-tola && r(ll,mm,nn,2) < y_bnd(kk,2)+tola && ...
                        r(ll,mm,nn,3) > z_bnd(kk,1)-tola && r(ll,mm,nn,3) < z_bnd(kk,2)+tola) 
                    boolean_tens(ll+(mm-1)*L+(nn-1)*L*M)=1;
                end
            end
        end
    end  
    idx_temp{kk}=zeros(sum(boolean_tens),3);
    idx_temp{kk}(:,1) = find(boolean_tens(:)==1);  
    idx_temp{kk}(:,2) = x_bnd(kk,3);  
    idx_temp{kk}(:,3) = x_bnd(kk,4); 
end
idx=cell2mat(idx_temp);  % 2nd dim indicates the conductor number
clear boolean_tens idx_temp x_bnd y_bnd z_bnd; 
disp(['Time for boolean loop ::: ',num2str(toc)])
% -------------------------------------------------------------------------
% Assign output
% -------------------------------------------------------------------------

epsilon_r = ones(L,M,N);
epsilon_r(idx(:,1)) = e_r;

sigma_e = zeros(L,M,N);
sigma_e(idx(:,1)) = s_e;

% -------------------------------------------------------------------------
% Assign the grid of conductors and visualize it 
% -------------------------------------------------------------------------

[L,M,N,~] = size(r);
grid_intcon = zeros(L,M,N,5);
grid_intcon(idx(:,1)) = r(idx(:,1));
grid_intcon(L*M*N+idx(:,1)) = r(L*M*N+idx(:,1));
grid_intcon(2*L*M*N+idx(:,1)) = r(2*L*M*N+idx(:,1));
grid_intcon(3*L*M*N+idx(:,1)) = idx(:,2);
grid_intcon(4*L*M*N+idx(:,1)) = idx(:,3);
% disp('Attention ::: Pad NaNs to grid_intcon for air voxels instead zeros! ')
% disp('This will avoid possible bug that appears when you have a voxel centered at (0,0,0)')

% plot geometry
%fl_plot_intcons=1;
if (fl_plot_intcons == 1)
    tic
    plot_boxes_of_grid(grid_intcon,dx);
    disp(['Time for visualizing the structure ::: ',num2str(toc)])
end


