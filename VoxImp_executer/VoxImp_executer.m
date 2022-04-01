clc; close all; clear all; format long e;
% -------------------------------------------------------------------------
%                  Add the Current Path to Workspace
% -------------------------------------------------------------------------
% parpool('local',31)
diary logfile.out
diary on
pre_define_the_path_for_folders
tot_time=tic;
% -------------------------------------------------------------------------
%                  Choose a geometry generator
% -------------------------------------------------------------------------
% Trans_line_weiping
% rec_spiral_blk_dia
trans_line_FastImp_Paper
% GND_x
% GND_TPH
% cross_buses_10x3
% cross_buses_10x6
% cross_buses_10x10
% test_small_bus
% parallel_buses_2x1
% L_shaped
% rec_spiral_thick
% Spiral_array
% -------------------------------------------------------------------------
%                  Flags and parameters
% -------------------------------------------------------------------------
% freq=[1e+08	 1.47368e+08 1.94737e+08 2.42105e+08 2.89474e+08 3.36842e+08 3.84211e+08 4.31579e+08 4.78947e+08 ...
%  5.26316e+08 5.73684e+08 6.21053e+08 6.68421e+08 7.15789e+08 7.63158e+08 8.10526e+08 8.57895e+08 9.05263e+08 9.52632e+08 1e+09];
% freq=[6 7.86206896551724 9.72413793103448 11.5862068965517 13.4482758620690 15.3103448275862 17.1724137931035 19.0344827586207 20.8965517241379 ...
%     22.7586206896552 24.6206896551724 26.4827586206897 28.3448275862069 30.2068965517241 32.0689655172414 33.9310344827586 35.7931034482759 ...
%     37.6551724137931 39.5172413793103 41.3793103448276 43.2413793103448 45.1034482758621 46.9655172413793 48.8275862068966 50.6896551724138  ...
%     52.5517241379310 54.4137931034483 56.2758620689655 58.1379310344828 60].*1e8;
% freq=[1e+09	1.08264e+09	1.1721e+09	1.26896e+09	1.37382e+09	1.48735e+09	1.61026e+09	1.74333e+09	1.88739e+09	2.04336e+09	...
% 2.21222e+09	2.39503e+09 2.59294e+09	2.80722e+09 3.0392e+09	3.29034e+09	3.56225e+09	3.85662e+09	4.17532e+09	4.52035e+09	...
% 4.8939e+09	5.29832e+09	5.73615e+09	6.21017e+09	6.72336e+09	7.27895e+09	7.88046e+09	8.53168e+09	9.23671e+09	1e+10];
% freq=[1 1.42857 1.85714 2.28571 2.71429 3.14286 3.57143 3.6 3.675 3.75 3.825 3.9 4 4.42857 4.85714 5.28571 5.71429 ...
%     6.14286 6.57143 7 7.42857 7.85714 8.28571 8.71429 9.14286 9.57143 10]*1e9;
freq=35e8;
% freq=[1e+07 1.26896e+07 1.61026e+07 2.04336e+07 2.59294e+07 3.29034e+07 4.17532e+07 5.29832e+07 6.72336e+07 8.53168e+07 ...
%     1.08264e+08 1.37382e+08 1.74333e+08 2.21222e+08 2.80722e+08 3.56225e+08 4.52035e+08 5.73615e+08 7.27895e+08 9.23671e+08 ...
%     1.1721e+09 1.48735e+09 1.88739e+09];
% freq=[5 5.51724137931035 6.03448275862069 6.55172413793104 7.06896551724138 7.58620689655172 8.10344827586207 8.62068965517242 ...
%     9.13793103448276 9.65517241379310 10.1724137931034 10.6896551724138 11.2068965517241 11.7241379310345 12.2413793103448 ...
%     12.7586206896552 13.2758620689655 13.7931034482759 14.3103448275862 14.8275862068966 15.3448275862069 15.8620689655172 ...
%     16.3793103448276 16.8965517241379 17.4137931034483 17.9310344827586 18.4482758620690 18.9655172413793 19.4827586206897 20].*1e9;
num_freq = length(freq);
epsa=[];% dielectric permittivity, for conductor case, set it to null []
epsb=1; % background permittivity
eps0=8.854187817e-12;
er = 0;  % epsilon_r of conductors
se=7.4e5; % conductivity of conductors 7.4e5 5.8e7;
% se=5.8e7;
inner_it = 50; outer_it = 20; tol=1e-6; % iterative solver inputs
fl_check_domain=0; % set to 1 for only plotting the structure (no simulation)
fl_check_geo=0; % set to 1 for only plotting the domain (no simulation)
fl_check_ports=0; % set to 1 for only plotting the port nodes (no simulation)
fl_check_domain_pnls=0; % 1 for plot the domain panels
fl_check_geo_pnls=0; % 1 for plot the bndry panels
plot_option=0; % see the options of plotting in Visualization part
freq_curr_plot=1e9; % frequency for plotting currents
fl_filling_circulant=0;
num_vox_in_blk=5; % size of block preconditioner
fl_Tuck=0;
tol_Tuck=1e-6;
% -------------------------------------------------------------------------
%                         Initialize stuff
% -------------------------------------------------------------------------

pre_print_out_inputs_generate_consts

% -------------------------------------------------------------------------
%                   Define domain and constitutive parameters
% -------------------------------------------------------------------------

% generate domain 3D grid
[r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);

% assign constitutive parameters
[idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams_Imp(r,Res,Cnt,Dims,Orients,er,se,fl_check_geo);

if (fl_check_domain == 1 || fl_check_geo == 1); return; end;

% -------------------------------------------------------------------------
%                 Define EM Vars/Constants and Domain Parameters
% -------------------------------------------------------------------------
pre_define_structure_params

% ------------------------------------------------------------------------
%                  Obtain Panel Coordinates and IDs
% -------------------------------------------------------------------------
tini_pre=tic;
disp('Generating panels of voxelized computational domain ...')
[comp_dom_panels]=new_computational_panels(dx,r);

disp('Generating panels of voxelized individual bodies ...')
% [geom_bndry_panels]=new_geom_panels(dx,grid_intcon,idx);
% [geom_bndry_panels]=port_geom_panels(dx,grid_intcon,idx,pnt_rght,pnt_lft,port_direction);
[geom_bndry_panels]=all_geom_panels(dx,grid_intcon,idx,num_vox_in_blk);

if fl_check_geo_pnls==1
    plot_panels(dx,geom_bndry_panels);
    return
end
if fl_check_domain_pnls==1
    plot_panels(dx,comp_dom_panels);
    return
end

geom_panels_cond=cell(num_conds,1);

for ii=1:num_conds
    geom_panels_cond{ii}=geom_bndry_panels(find(geom_bndry_panels(:,6)==ii),:);
    geom_panels_cond{ii}(:,8)=Eps_inout(ii,1);
    geom_panels_cond{ii}(:,9)=Eps_inout(ii,2);
end
if num_diels>0
    geom_panels_diel=cell(num_diels,1);
    for ii=1:num_diels
        geom_panels_diel{ii}=geom_bndry_panels(find(geom_bndry_panels(:,7)==ii),:);
        geom_panels_diel{ii}(:,8)=Eps_inout(ii,1);
        geom_panels_diel{ii}(:,9)=Eps_inout(ii,2);
    end
else
    geom_panels_diel=[];
end

fl_plot_pnl_orients=0;
geom_bndry_panels=[geom_panels_cond; geom_panels_diel];
geom_bndry_panels=mingyu_pre_correct_pnl_orients(num_conds,num_diels,geom_bndry_panels,Cnt,dx,fl_plot_pnl_orients);

% defining 1 or 2 for activating structures for conductors or conductors w/ diel
% We have one data structure that has combined conductor info and one data
% structure that wil have combined dielectric info.
if (num_diels == 0)
    fl_st_diel=1;
else
    fl_st_diel=num_diels+1;
end

num_panels=cell(fl_st_diel,1);
ids_panels=cell(fl_st_diel,1);
num_tot_panels=0;
geom_bndry_panels_cond=cell2mat(geom_panels_cond);
geom_bndry_panels_cond=sortrows(geom_bndry_panels_cond,5);
for kk=1:1
    [num_panels{kk},ids_panels{kk}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels_cond);
    num_tot_panels=num_tot_panels+num_panels{kk};
end
if num_diels>0
    for kk=1:num_diels
        [num_panels{kk+1},ids_panels{kk+1}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_panels_diel{kk});
        num_tot_panels=num_tot_panels+num_panels{kk+1};
    end
end
clear comp_dom_panels geom_panels_cond geom_panels_diel

if (num_diels == 0)
    ids_panels{2}=cell(6,1);
end

geom_bndry_panels2=cell(fl_st_diel,1);
geom_bndry_panels2{1}=geom_bndry_panels_cond;
clear geom_bndry_panels_cond

if(num_diels > 0)% just add at the end
    for ii=1:num_diels
        geom_bndry_panels2{ii+1} = geom_bndry_panels{num_conds+ii};
    end
end
clear  geom_bndry_panels

geom_bndry_panels=zeros(num_tot_panels,11);
st_ind=1; end_ind=0;
for kk=1:fl_st_diel
    end_ind=end_ind+num_panels{kk};
    geom_bndry_panels(st_ind:end_ind,1:11)=geom_bndry_panels2{kk}(:,1:11);
    st_ind=end_ind+1;
end
clear  geom_bndry_panels2

diel_pnl_id_locs=cell(num_diels,1);diel_pnl_ids=cell(num_diels,1);
for ii=1:num_diels
    diel_pnl_id_locs{ii}=cell(6,1); diel_pnl_ids{ii}=cell(6,1);
    if (num_diels > 0) % currently num_diels could be max 1
        % information for +/- direction directed panels from combined ids_panels
        diel_pnl_id_locs{ii}{1}=find(ids_panels{1+ii}{4}(:)>0);
        diel_pnl_ids{ii}{1}=ids_panels{1+ii}{1}(diel_pnl_id_locs{ii}{1});
        diel_pnl_id_locs{ii}{4}=find(ids_panels{1+ii}{4}(:)<0);
        diel_pnl_ids{ii}{4}=ids_panels{1+ii}{1}(diel_pnl_id_locs{ii}{4});
        
        diel_pnl_id_locs{ii}{2}=find(ids_panels{1+ii}{5}(:)>0);
        diel_pnl_ids{ii}{2}=ids_panels{1+ii}{2}(diel_pnl_id_locs{ii}{2});
        diel_pnl_id_locs{ii}{5}=find(ids_panels{1+ii}{5}(:)<0);
        diel_pnl_ids{ii}{5}=ids_panels{1+ii}{2}(diel_pnl_id_locs{ii}{5});
        
        diel_pnl_id_locs{ii}{3}=find(ids_panels{1+ii}{6}(:)>0);
        diel_pnl_ids{ii}{3}=ids_panels{1+ii}{3}(diel_pnl_id_locs{ii}{3});
        diel_pnl_id_locs{ii}{6}=find(ids_panels{1+ii}{6}(:)<0);
        diel_pnl_ids{ii}{6}=ids_panels{1+ii}{3}(diel_pnl_id_locs{ii}{6});
        
        % diel_pnl_id_locs stores the locations of elements of diel_pnl_ids which
        % are positive/negative x, y, and z directed panels.
        % Its entries {1} +x-directed, {2} +y-directed, {3} +z-directed, {4}
        % -x-directed, {5} -y directed, and {6} -z-directed
    end
end
% Attention: the following lines were added to distinguish x-, y-, z- aligned
% conductor and dielectric panels

num_x_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 1)); num_y_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 2)); num_z_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 3));
num_x_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0));
num_y_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0));
num_z_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0));

if num_diels>0
    num_x_aligned_diel=zeros(num_diels,1);num_y_aligned_diel=zeros(num_diels,1);num_z_aligned_diel=zeros(num_diels,1);
    for ii=1:num_diels
        num_x_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == ii));
        num_y_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == ii));
        num_z_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == ii));
    end
end

inds_glob=zeros(num_diels+1,2);
inds_glob(1,:)=[1 num_x_aligned_cond];
inds_glob(num_diels+2,:) = [num_x_aligned_all+1 num_x_aligned_all+num_y_aligned_cond];
inds_glob(num_diels*2+3,:) = [num_x_aligned_all+num_y_aligned_all+1 num_x_aligned_all+num_y_aligned_all+num_z_aligned_cond]; % conductor part

if num_diels>0
    for ii=1:num_diels
        inds_glob(1+ii,:)=[inds_glob(ii,2)+1 inds_glob(ii,2)+sum(num_x_aligned_diel(ii))];
        inds_glob(num_diels+2+ii,:)=[inds_glob(num_diels+2+ii-1,2)+1 inds_glob(num_diels+2+ii-1,2)+sum(num_y_aligned_diel(ii))];
        inds_glob(num_diels*2+3+ii,:)=[inds_glob(num_diels*2+3+ii-1,2)+1 inds_glob(num_diels*2+3+ii-1,2)+sum(num_z_aligned_diel(ii))];
    end
end
tend = toc(tini_pre);
sim_preproc=tend;
disp(['Total time for generate panels ::: ' ,num2str(sim_preproc)]);
s=whos('geom_bndry_panels');
memory_bndry_panels=s.bytes/1024^2;
memory_stage_preprocessing=memory_bndry_panels
% index_bndry_node=geom_bndry_panels(:,10);
num_diel=0;
[blk_pre2,ids_pre,ord_pre]=implement_preconditioner_ver2nn(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel,freq*dx^4);

% ------------------------------------------------------------------------
%                  Obtain Nodal Incidence Matrix
% -------------------------------------------------------------------------
tini = tic;
% [Ae_original,nodeid_lft,nodeid_rght,nodeid_wlcond,Ae_only_leaving,Ae_only_entering_bndry] = lse_compute_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft,pnt_rght,pnt_well_cond,fl_check_ports);
[AL,nodeid_lft,nodeid_rght,nodeid_wlcond,AL_only_leaving,AL_only_entering_bndry,flag_node_cond_no] = compute_AL_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft, pnt_rght,pnt_well_cond,fl_check_ports);

if (fl_check_ports == 1); return; end;
tend = toc(tini);
disp(['Time for generating Ae mat & finding IDs of port nodes::: ' ,num2str(tend)]);

clear grid_intcon idx;
% ------------------------------------------------------------------------
%              RHS & Circulant Tensors
% -------------------------------------------------------------------------
disp('-----------------------------------------------------')
disp(['Precomputing LSE data structures...'])
N_c=size(AL,2); Nn=size(AL,1); tot_num=N_c+2*Nn;
rhs_vect=zeros(Nn,num_ports); % I is excitation
RHS=zeros(tot_num,num_ports);
for port_no=1:num_ports
    % ------------------------------------------------------------------------
    %              Assign Excitation and Ground Nodes
    % ------------------------------------------------------------------------
    [nodeid_4_grnd,nodeid_4_injectcurr,nodeid_4_grnd_port]=lse_assign_exc_grnd_nodes(nodeid_lft,nodeid_rght,nodeid_wlcond,num_ports,port_no);
    % ------------------------------------------------------------------------
    %                  Set Excitation (V)
    % ------------------------------------------------------------------------
    rhs_vect(nodeid_4_injectcurr,port_no)=-1/(wid/dx*hei/dx);
    rhs_vect(nodeid_4_grnd_port,port_no)=1/(wid/dx*hei/dx);
    [RHS(N_c+Nn+1:N_c+2*Nn,port_no)] = rhs_vect(:,port_no); % I is source
   
        index=1:Nn;
        ones_vect=ones(length(index),1);
        Aq=sparse(index,index,ones_vect,Nn,Nn);
   
    A_ind=[AL Aq];
    
    
    % ------------------------------------------------------------------------
    %         Generate Circulant Tensors
    % -------------------------------------------------------------------------
    
    if (port_no == 1)
        tini = tic;
        if fl_filling_circulant==0
            load('prestore_2000_10_30.mat');% according to the structure size to choose prestored Toeplitz
%             load('prestore_250_250_100.mat')
%             load('prestore_3000_2000_10.mat')
            %             load('prestore_500_500_50.mat');% according to the structure size to choose prestored Toeplitz
%             load('prestore_8000_20_60.mat');% according to the structure size to choose prestored Toeplitz
            [fN_ch_all]=retrieval_circulant_T_cap(dx,L,M,N,Eps_inout,tol_Tuck,fl_Tuck,T_cap);
            cap_dia1=1/fN_ch_all{1,1,1}(1,1,1);
            fl_no_fft=0;
            [fN_all,st_sparse_precon] = retrieval_circulant_T_henry(dx,1,L,M,N,fl_no_fft,T_henry,fl_Tuck,tol_Tuck);
        else
            [fN_ch_all]=generate_circulant_tensor_charge(dx,L,M,N,Eps_inout,1e-4,0);
            if (num_freq == 1)
                fl_no_fft=0;
                [fN_all,st_sparse_precon] = lse_generate_circulant_tensor_henry(dx,ko,L,M,N,fl_no_fft);
            else
                fl_no_fft=1;
                [fN_all2,st_sparse_precon2] = lse_generate_circulant_tensor_henry(dx,1,L,M,N,fl_no_fft);
                %                 load fN.mat
                % note multiply fN_all and st_sparse_precon with ko^2 and compute its FFT
            end
        end
        tend = toc(tini);
        disp(['Total time for getting circulant tensor ::: ' ,num2str(tend)]);
        
    end
end

disp(['Done... Precomputing LSE data structures'])
disp('-----------------------------------------------------')
% ------------------------------------------------------------------------
%              Slove LSE
% -------------------------------------------------------------------------
disp(['Solving LSEs ...'])
Y_mat=zeros(num_ports,num_ports,num_freq);
Y_mat2=zeros(num_ports,num_ports,num_freq);
Z_mat=zeros(num_ports,num_ports,num_freq);
R_jL_mat=zeros(num_ports,num_ports,num_freq);

for freq_no=1:num_freq
    tinisim = tic;
    if (num_freq > 1)
        freq = freq_all(freq_no);
    end
    EMconstants
    disp('-----------------------------------------------------')
    disp(['Simulation for frequency : ',num2str(freq),' started! ', 'freq pnt: ',num2str(freq_no), ' / ', num2str(num_freq)])
    % setting new constitutive parameters for new freq
    Mr = epsilon_r - 1j*sigma_e/(eo*omega); % permittivity
    Mc = Mr - 1.0; % susceptibility
    OneoverMc = 1.0 ./ Mc; % one over susceptibility
    
    for port_no=1:1%num_ports
        disp(['Solving for port # ',num2str(port_no), ' ...'])
        % ------------------------------------------------------------------------
        %     Solve Linear System of Equations Iteratively
        % -------------------------------------------------------------------------
        dim_ind=[N_c Nn];
        if (port_no == 1)
            % prepare the preconditioner
            tinisim = tic;
            blk_pre=cell(size(ids_pre,1),1);
            for ll=1:size(ids_pre,1)
                blk_pre{ll}=blk_pre2{ll}*(2*pi*freq*dx^4);
            end
            cap_dia=cap_dia1*(2*pi*freq*dx^4);
            lse_sparse_precon_prepare(dx,freq,OneoverMc,idxS3,st_sparse_precon*ko^2,AL,dim_ind,blk_pre,ids_pre);
            disp(['Time for prepare precond ::: ' ,num2str(toc(tinisim))]);
        end
        % Define the handle for matvect

        fACPU   = @(J)lse_matvect_mult1(J, fN_all,fN_ch_all, A_ind,dim_ind, OneoverMc, dx, freq, idxS5,1,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,num_diels,ko,blk_pre,ids_pre,ord_pre,fl_Tuck);
        tini = tic;
        disp(['Iterative solution started ... '])
        [rhs_vect_sparse_precon]=lse_sparse_precon_multiply(RHS(:,port_no),AL,Aq,blk_pre,ids_pre,ord_pre);
%         [rhs_vect_sparse_precon]=lse_sparse_precon_multiply_dia(RHS(:,port_no),AL,Aq);
        [x, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J), rhs_vect_sparse_precon, inner_it, tol, outer_it);
        tend = toc(tini);
        disp(['Total time for iterative solution ::: ' ,num2str(tend)]);
        disp(['Done... Iterative solution'])
        
        if (abs(freq_curr_plot-freq)<1e-12 && port_no == 1)
            x_backup = x;
        end
        % ------------------------------------------------------------------------
        %     Compute the Currents on Port Nodes and the Column in Ymat
        % -------------------------------------------------------------------------
        currs_port_yparams=zeros(num_ports,1);
        for kk=1:num_ports
            currs_port_yparams(kk,1)=sum(RHS(:,kk).*x);
        end
        
        Y_mat(:,port_no,freq_no)=currs_port_yparams(:,1);
        disp(['Done... Solving for port # ',num2str(port_no)])
        
    end
    disp('-----------------------------------------------------')
    disp(['Done... Simulation for frequency : ',num2str(freq),' freq pnt: ',num2str(freq_no), ' / ', num2str(num_freq)])
    Z_mat(:,:,freq_no)=(squeeze(Y_mat(:,:,freq_no))); % I is source
    R_jL_mat(:,:,freq_no)=(real(squeeze(Z_mat(:,:,freq_no))))+sqrt(-1)*(imag(Z_mat(:,:,freq_no)));%/(2*pi*freq));
    
end
RJLMAT=squeeze(R_jL_mat);
ZZ=sqrt(real(RJLMAT).^2+imag(RJLMAT).^2);

toc(tot_time)


diary off
