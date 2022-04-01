% clc; close all; clear all; format long e;
% % -------------------------------------------------------------------------
% %                  Add the Current Path to Workspace
% % -------------------------------------------------------------------------
% % parpool('local',31)
% pre_define_the_path_for_folders
% 
% % -------------------------------------------------------------------------
% %                  Inputs for Simulation
% % -------------------------------------------------------------------------
% % freq=[1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]*1e9;
% freq=[1e+09	1.08264e+09	1.1721e+09	1.26896e+09	1.37382e+09	1.48735e+09	1.61026e+09	1.74333e+09	1.88739e+09	2.04336e+09	...
% 2.21222e+09	2.39503e+09 2.59294e+09	2.80722e+09 3.0392e+09	3.29034e+09	3.56225e+09	3.85662e+09	4.17532e+09	4.52035e+09	...
% 4.8939e+09	5.29832e+09	5.73615e+09	6.21017e+09	6.72336e+09	7.27895e+09	7.88046e+09	8.53168e+09	9.23671e+09	1e+10];
% num_freq = length(freq);
% epsa=[];% dielectric permittivity, for conductor case, set it to null []
% epsb=1; % background permittivity
% eps0=8.854187817e-12;
% er = 0;  % epsilon_r of conductors
% se=7.4e5; % conductivity of conductors 7.4e5 5.8e7;
% inner_it = 150; outer_it = 10; tol=1e-6; % iterative solver inputs
Res = 10e-6; % voxel size (deltax)
% fl_check_domain=0; % set to 1 for only plotting the structure (no simulation)
% fl_check_geo=0; % set to 1 for only plotting the domain (no simulation)
% fl_check_ports=0; % set to 1 for only plotting the port nodes (no simulation)
% fl_check_domain_pnls=0; % 1 for plot the domain panels
% fl_check_geo_pnls=0; % 1 for plot the bndry panels
% plot_option=0; % see the options of plotting in Visualization part
% freq_curr_plot=1e9; % frequency for plotting currents
% fl_filling_circulant=0;
% num_vox_in_blk=5;
% fl_Vsour=2; % 1 for V source, 2 for I source
% -------------------------------------------------------------------------
%                  Inputs for the Structure
% -------------------------------------------------------------------------
% We only need centers (Cnt), dimensions (Dims), and orientations (Orients)
% of the conductors at the end of this part.

% inputs for generating conductors with specified lengths and widths of arms
num_conds = 1; % number of conductors
num_diels = 0;
num_ports = 1; % number of ports

len1=20000e-6; wid1=50e-6; hei1=50e-6;
len2=50e-6; wid2=10e-6; hei2=50e-6;

cen_temp1=[len1/2 wid1/2 hei1/2 1 0]; %4th col:cond numbering, 0 for GND or diel, %5th col: diel numbering, 0 for no diel
cen_temp2=[len1/2 wid1/2 hei1+len2+hei1/2 1 0];
cen_temp3=[len1-wid2/2 wid1/2 hei1+len2/2 1 0];


Cnt = [cen_temp1 ;cen_temp2;cen_temp3;]; % centers of conductors
Dims_tmp1 = [len1 wid1 hei1]; Dims_tmp2 = [len1 wid1 hei1]; Dims_tmp3 = [len2 hei2 wid2 ]; % dimensions of conductors(L(x),W(y),H(z))
Orients_tmp=['x';'x';'z';]; % orientations of conductors %
Dims=[Dims_tmp1;Dims_tmp2;Dims_tmp3;];
Orients=[Orients_tmp];

Eps_inout=zeros(num_conds+num_diels,2);% inner and outer eps
for ii=1:num_conds
    Eps_inout(ii,:)=[0 1];
end
% -------------------------------------------------------------------------
%                  Input for Computational Domain
% -------------------------------------------------------------------------
% At the end of this part, we only need bbox_min(3) and bbox_max(3) vectors
% define computational domain or bounding box enclosing the structure
bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[len1 wid1 2*hei1+len2]+1e-12; % max coordinates of bbox

% -------------------------------------------------------------------------
%                  Input for Ports
% -------------------------------------------------------------------------
% At the end of this part, we need structures pnt_lft{xx} and pnt_rght{xx} which
% contains the coordinates of nodes on both sides of xxth port
port_direction=[1 1];%normal direction of left and right port. 1-x 2-y 3-z
% defining the nodes in first port
pnt_lft=cell(num_ports,1); %exc
pnt_rght=cell(num_ports,1); %gnd
pnt_lft{1}=zeros(round(wid1/Res)*round(hei1/Res),4);
pnt_rght{1}=zeros(round(wid1/Res)*round(hei1/Res),4);
dum=1;
% for kk=1:round(wid/Res)
%     for ll=1:round(hei/Res)
%         pnt_rght{1}(dum,1:3)=[len2+len15-wid (2*kk-1)*(0.5*Res)+2*wid+2*dd (2*ll-1)*(0.5*Res)]; % points on which excitation defined
%         pnt_lft{1}(dum,1:3)=[len2+len15-wid (2*kk-1)*(0.5*Res)+len3-len1 (2*ll-1)*(0.5*Res)]; % points on which ground defined
%         dum=dum+1;
%     end
% end
for kk=1:round(wid1/Res)
    for ll=1:round(hei1/Res)
        pnt_rght{1}(dum,1:4)=[0 (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res) 1]; 
        pnt_lft{1}(dum,1:4)=[0  (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res)+hei1+len2 1]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
% defining nodes connected ground if conductors without ports exist; if
% there is no, then leave as a empty array.

%pnt_well_cond=[pnt_lft{1}(1,1) pnt_lft{1}(1,2) + dist_btw_conds pnt_lft{1}(1,3);];
pnt_well_cond=[];
wid=wid1; hei=hei1;
% % -------------------------------------------------------------------------
% %                         Initialize stuff
% % -------------------------------------------------------------------------
% 
% pre_print_out_inputs_generate_consts
% 
% % -------------------------------------------------------------------------
% %                   Define domain and constitutive parameters
% % -------------------------------------------------------------------------
% 
% % generate domain 3D grid
% [r] = generategridfrombbox(Res,[bbox_min(1) bbox_max(1)],[bbox_min(2) bbox_max(2)],[bbox_min(3) bbox_max(3)],fl_check_domain);
% 
% % assign constitutive parameters
% [idx,epsilon_r,sigma_e,grid_intcon] = intcon_constparams_Imp(r,Res,Cnt,Dims,Orients,er,se,fl_check_geo);
% 
% if (fl_check_domain == 1 || fl_check_geo == 1); return; end;
% 
% % -------------------------------------------------------------------------
% %                 Define EM Vars/Constants and Domain Parameters
% % -------------------------------------------------------------------------
% pre_define_structure_params
% 
% % % ------------------------------------------------------------------------
% % %                  Obtain Nodal Incidence Matrix
% % % -------------------------------------------------------------------------
% % tini = tic;
% % % [Ae_original,nodeid_lft,nodeid_rght,nodeid_wlcond,Ae_only_leaving,Ae_only_entering_bndry] = lse_compute_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft,pnt_rght,pnt_well_cond,fl_check_ports);
% % port_direction=[1 1];%normal direction of left and right port. 1-x 2-y 3-z
% % [Ae,Ai_original1,nodeid_lft,nodeid_rght,nodeid_wlcond,Ae_only_leaving,Ae_only_entering_bndry] = compute_Ai_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft,pnt_rght,port_direction,fl_check_ports);
% % Ai_original1(all(Ai_original1==0,2),:)=[];
% % if (fl_check_ports == 1); return; end;
% % tend = toc(tini);
% % disp(['Time for generating Ae mat & finding IDs of port nodes::: ' ,num2str(tend)]);
% 
% 
% % ------------------------------------------------------------------------
% %                  Obtain Panel Coordinates and IDs
% % -------------------------------------------------------------------------
% tini_pre=tic;
% disp('Generating panels of voxelized computational domain ...')
% [comp_dom_panels]=new_computational_panels(dx,r);
% clear r
% 
% disp('Generating panels of voxelized individual bodies ...')
% % [geom_bndry_panels]=new_geom_panels(dx,grid_intcon,idx);
% % [geom_bndry_panels]=port_geom_panels(dx,grid_intcon,idx,pnt_rght,pnt_lft,port_direction);
% [geom_bndry_panels]=all_geom_panels(dx,grid_intcon,idx,num_vox_in_blk);
% 
% if fl_check_geo_pnls==1
%     plot_panels(dx,geom_bndry_panels);
%     return
% end
% if fl_check_domain_pnls==1
%     plot_panels(dx,comp_dom_panels);
%     return
% end
% 
% geom_panels_cond=cell(num_conds,1);
% 
% for ii=1:num_conds
%     geom_panels_cond{ii}=geom_bndry_panels(find(geom_bndry_panels(:,6)==ii),:);
%     geom_panels_cond{ii}(:,8)=Eps_inout(ii,1);
%     geom_panels_cond{ii}(:,9)=Eps_inout(ii,2);
% end
% if num_diels>0
%     geom_panels_diel=cell(num_diels,1);
%     for ii=1:num_diels
%         geom_panels_diel{ii}=geom_bndry_panels(find(geom_bndry_panels(:,7)==ii),:);
%         geom_panels_diel{ii}(:,8)=Eps_inout(ii,1);
%         geom_panels_diel{ii}(:,9)=Eps_inout(ii,2);
%     end
% else
%     geom_panels_diel=[];
% end
% 
% fl_plot_pnl_orients=0;
% geom_bndry_panels=[geom_panels_cond; geom_panels_diel];
% geom_bndry_panels=mingyu_pre_correct_pnl_orients(num_conds,num_diels,geom_bndry_panels,Cnt,dx,fl_plot_pnl_orients);
% 
% % defining 1 or 2 for activating structures for conductors or conductors w/ diel
% % We have one data structure that has combined conductor info and one data
% % structure that wil have combined dielectric info.
% if (num_diels == 0)
%     fl_st_diel=1;
% else
%     fl_st_diel=num_diels+1;
% end
% 
% num_panels=cell(fl_st_diel,1);
% ids_panels=cell(fl_st_diel,1);
% num_tot_panels=0;
% geom_bndry_panels_cond=cell2mat(geom_panels_cond);
% geom_bndry_panels_cond=sortrows(geom_bndry_panels_cond,5);
% for kk=1:1
%     [num_panels{kk},ids_panels{kk}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_bndry_panels_cond);
%     num_tot_panels=num_tot_panels+num_panels{kk};
% end
% if num_diels>0
%     for kk=1:num_diels
%         [num_panels{kk+1},ids_panels{kk+1}]=pre_conn_domain_and_bndry_panels(comp_dom_panels,geom_panels_diel{kk});
%         num_tot_panels=num_tot_panels+num_panels{kk+1};
%     end
% end
% clear comp_dom_panels geom_panels_cond geom_panels_diel
% 
% if (num_diels == 0)
%     ids_panels{2}=cell(6,1);
% end
% 
% geom_bndry_panels2=cell(fl_st_diel,1);
% geom_bndry_panels2{1}=geom_bndry_panels_cond;
% clear geom_bndry_panels_cond
% 
% if(num_diels > 0)% just add at the end
%     for ii=1:num_diels
%         geom_bndry_panels2{ii+1} = geom_bndry_panels{num_conds+ii};
%     end
% end
% clear  geom_bndry_panels
% 
% geom_bndry_panels=zeros(num_tot_panels,11);
% st_ind=1; end_ind=0;
% for kk=1:fl_st_diel
%     end_ind=end_ind+num_panels{kk};
%     geom_bndry_panels(st_ind:end_ind,1:11)=geom_bndry_panels2{kk}(:,1:11);
%     st_ind=end_ind+1;
% end
% clear  geom_bndry_panels2
% 
% diel_pnl_id_locs=cell(num_diels,1);diel_pnl_ids=cell(num_diels,1);
% for ii=1:num_diels
%     diel_pnl_id_locs{ii}=cell(6,1); diel_pnl_ids{ii}=cell(6,1);
%     if (num_diels > 0) % currently num_diels could be max 1
%         % information for +/- direction directed panels from combined ids_panels
%         diel_pnl_id_locs{ii}{1}=find(ids_panels{1+ii}{4}(:)>0);
%         diel_pnl_ids{ii}{1}=ids_panels{1+ii}{1}(diel_pnl_id_locs{ii}{1});
%         diel_pnl_id_locs{ii}{4}=find(ids_panels{1+ii}{4}(:)<0);
%         diel_pnl_ids{ii}{4}=ids_panels{1+ii}{1}(diel_pnl_id_locs{ii}{4});
%         
%         diel_pnl_id_locs{ii}{2}=find(ids_panels{1+ii}{5}(:)>0);
%         diel_pnl_ids{ii}{2}=ids_panels{1+ii}{2}(diel_pnl_id_locs{ii}{2});
%         diel_pnl_id_locs{ii}{5}=find(ids_panels{1+ii}{5}(:)<0);
%         diel_pnl_ids{ii}{5}=ids_panels{1+ii}{2}(diel_pnl_id_locs{ii}{5});
%         
%         diel_pnl_id_locs{ii}{3}=find(ids_panels{1+ii}{6}(:)>0);
%         diel_pnl_ids{ii}{3}=ids_panels{1+ii}{3}(diel_pnl_id_locs{ii}{3});
%         diel_pnl_id_locs{ii}{6}=find(ids_panels{1+ii}{6}(:)<0);
%         diel_pnl_ids{ii}{6}=ids_panels{1+ii}{3}(diel_pnl_id_locs{ii}{6});
%         
%         % diel_pnl_id_locs stores the locations of elements of diel_pnl_ids which
%         % are positive/negative x, y, and z directed panels.
%         % Its entries {1} +x-directed, {2} +y-directed, {3} +z-directed, {4}
%         % -x-directed, {5} -y directed, and {6} -z-directed
%     end
% end
% % Attention: the following lines were added to distinguish x-, y-, z- aligned
% % conductor and dielectric panels
% 
% num_x_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 1)); num_y_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 2)); num_z_aligned_all=length(find(abs(geom_bndry_panels(:,4)) == 3));
% num_x_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == 0));
% num_y_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == 0));
% num_z_aligned_cond=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == 0));
% 
% if num_diels>0
%     num_x_aligned_diel=zeros(num_diels,1);num_y_aligned_diel=zeros(num_diels,1);num_z_aligned_diel=zeros(num_diels,1);
%     for ii=1:num_diels
%         num_x_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 1 & geom_bndry_panels(:,7) == ii));
%         num_y_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 2 & geom_bndry_panels(:,7) == ii));
%         num_z_aligned_diel(ii)=length(find(abs(geom_bndry_panels(:,4)) == 3 & geom_bndry_panels(:,7) == ii));
%     end
% end
% 
% inds_glob=zeros(num_diels+1,2);
% inds_glob(1,:)=[1 num_x_aligned_cond];
% inds_glob(num_diels+2,:) = [num_x_aligned_all+1 num_x_aligned_all+num_y_aligned_cond];
% inds_glob(num_diels*2+3,:) = [num_x_aligned_all+num_y_aligned_all+1 num_x_aligned_all+num_y_aligned_all+num_z_aligned_cond]; % conductor part
% 
% if num_diels>0
%     for ii=1:num_diels
%         inds_glob(1+ii,:)=[inds_glob(ii,2)+1 inds_glob(ii,2)+sum(num_x_aligned_diel(ii))];
%         inds_glob(num_diels+2+ii,:)=[inds_glob(num_diels+2+ii-1,2)+1 inds_glob(num_diels+2+ii-1,2)+sum(num_y_aligned_diel(ii))];
%         inds_glob(num_diels*2+3+ii,:)=[inds_glob(num_diels*2+3+ii-1,2)+1 inds_glob(num_diels*2+3+ii-1,2)+sum(num_z_aligned_diel(ii))];
%     end
% end
% tend = toc(tini_pre);
% sim_preproc=tend;
% disp(['Total time for generate panels ::: ' ,num2str(sim_preproc)]);
% s=whos('geom_bndry_panels');
% memory_bndry_panels=s.bytes/1024^2;
% memory_stage_preprocessing=memory_bndry_panels
% % index_bndry_node=geom_bndry_panels(:,10);
% num_diel=0;
% [blk_pre2,ids_pre,ord_pre]=implement_preconditioner_ver2nn(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel,freq*dx^4);
% 
% % ------------------------------------------------------------------------
% %                  Obtain Nodal Incidence Matrix
% % -------------------------------------------------------------------------
% tini = tic;
% % [Ae_original,nodeid_lft,nodeid_rght,nodeid_wlcond,Ae_only_leaving,Ae_only_entering_bndry] = lse_compute_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft,pnt_rght,pnt_well_cond,fl_check_ports);
% [AL,nodeid_lft,nodeid_rght,nodeid_wlcond,AL_only_leaving,AL_only_entering_bndry] = compute_AL_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft, pnt_rght,port_direction,fl_check_ports);
% 
% % [Ae,Ai_original1,nodeid_lft,nodeid_rght,nodeid_wlcond,Ae_only_leaving,Ae_only_entering_bndry] = compute_Ai_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft,pnt_rght,port_direction,index_bndry_node,fl_check_ports);
% % Ai_original1(all(Ai_original1==0,2),:)=[];
% if (fl_check_ports == 1); return; end;
% tend = toc(tini);
% disp(['Time for generating Ae mat & finding IDs of port nodes::: ' ,num2str(tend)]);
% 
% clear grid_intcon idx;
% % ------------------------------------------------------------------------
% %              RHS & Circulant Tensors
% % -------------------------------------------------------------------------
% disp('-----------------------------------------------------')
% disp(['Precomputing LSE data structures...'])
% % N_c=size(Ae,2); N_in=size(Ai_original1,1); N_out=size(Ae,1); tot_num=N_c+2*N_out+N_in;
% % temp=size([Ae;Ai_original1],1)+size([Ae;Ai_original1],2);
% % % rhs_vect=zeros(temp,num_ports); % V is excitation
% % rhs_vect=zeros(N_out,num_ports); % I is excitation
% % RHS=zeros(tot_num,num_ports);
% 
% N_c=size(AL,2); Nn=size(AL,1); tot_num=N_c+2*Nn;
% rhs_vect=zeros(Nn,num_ports); % I is excitation
% RHS=zeros(tot_num,num_ports);
% for port_no=1:num_ports
%     % ------------------------------------------------------------------------
%     %              Assign Excitation and Ground Nodes
%     % ------------------------------------------------------------------------
%     
%     [nodeid_4_grnd,nodeid_4_injectcurr]=lse_assign_exc_grnd_nodes(nodeid_lft,nodeid_rght,nodeid_wlcond,num_ports,port_no);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Aind
%     %     %     index=(find(geom_bndry_panels(:,6)==port_no));
%     % %     ones_vect=ones(length(index),1);
%     %     index1=(nodeid_4_injectcurr);
%     %     index2=(nodeid_4_grnd);
%     %     index=[index1; index2];
%     %     ones_vect1=ones(length(index1),1);
%     %     ones_vect2=ones(length(index2),1);
%     %     ones_vect=[ones_vect1; ones_vect2];
%     %     Aq=sparse(index,index,ones_vect,N_out,N_out);
%     %     Zeros_ind=sparse(N_in,N_out);
%     %     A_ind=[Ae Aq; Ai_original1 Zeros_ind];
%     %
%     %     Zeros_q=sparse(N_out,N_out);
%     %     A_ind1=[Ae Zeros_q; Ai_original1 Zeros_ind];
%     
%     sc=size(AL,2)/5;
%     %     index1=(nodeid_4_injectcurr);
%     %     index2=(nodeid_4_grnd);
%     %     index=[index1; index2];
%     %     ones_vect1=ones(length(index1),1);
%     %     ones_vect2=ones(length(index2),1);
%     %     ones_vect=[ones_vect1; ones_vect2];
%     index=1:Nn;
%     %     index=find(geom_bndry_panels(:,11)==1);
%     ones_vect=ones(length(index),1);
%     
%     Aq=sparse(index,index,ones_vect,Nn,Nn);
%     %     Zeros_ind=sparse(Nn,Nn);
%     A_ind=[AL Aq];
%     %     Ai = Ai_original1;
%     
%     % ------------------------------------------------------------------------
%     %                  Set Excitation (V)
%     % ------------------------------------------------------------------------
%     if fl_Vsour==2
%         rhs_vect(nodeid_4_injectcurr,port_no)=-1/(wid1/dx*hei1/dx);
%         rhs_vect(nodeid_4_grnd,port_no)=1/(wid1/dx*hei1/dx);
%         [RHS(N_c+Nn+1:N_c+2*Nn,port_no)] = rhs_vect(:,port_no); % I is source
%     else
%         [RHS(:,port_no)] = lse_compute_rhs_vector(A_ind,nodeid_4_injectcurr); % V is source
%     end
%     
%     % ------------------------------------------------------------------------
%     %         Generate Circulant Tensors
%     % -------------------------------------------------------------------------
%     
%     if (port_no == 1)
%         tini = tic;
%         if fl_filling_circulant==0
% %             size_up_lim=max([L,M,N]);
% %             if size_up_lim <=50
% %                 load('fN_prestored_data_50.mat'); % according to the structure size to choose prestored Toeplitz
% %             elseif size_up_lim>50 && size_up_lim <=100
% %                 load('fN_prestored_data_100.mat'); % according to the structure size to choose prestored Toeplitz
% %             elseif size_up_lim>100 && size_up_lim <=200
% %                 load('fN_prestored_data_200.mat'); % according to the structure size to choose prestored Toeplitz
% %             elseif size_up_lim>200 && size_up_lim <=300
% %                 load('fN_prestored_data_300.mat'); % according to the structure size to choose prestored Toeplitz
% %             elseif size_up_lim>300 && size_up_lim <=400
% %                 load('fN_prestored_data_400.mat'); % according to the structure size to choose prestored Toeplitz
% %             else
% %                 load('fN_prestored_data_500.mat'); % according to the structure size to choose prestored Toeplitz
% %             end
%             load prestore_Toep_trans.mat
%             [fN_ch_all]=retrival_circulant_cap(dx,L,M,N,Eps_inout,1e-6,1,trans_ch);
%             
%             if (num_freq == 1)
%                 fl_no_fft=0;
%                 [fN_all,st_sparse_precon] = retrieval_circulant_henry(dx,ko,L,M,N,fl_no_fft,trans_hen);
%             else
%                 fl_no_fft=1;
%                 [fN_all2,st_sparse_precon2] = retrieval_circulant_henry(dx,1,L,M,N,fl_no_fft,trans_hen);
%                 % note multiply fN_all and st_sparse_precon with ko^2 and compute its FFT
%             end
%         else
%             [fN_ch_all]=generate_circulant_tensor_charge(dx,L,M,N,Eps_inout,1e-4,0);
%             if (num_freq == 1)
%                 fl_no_fft=0;
%                 [fN_all,st_sparse_precon] = lse_generate_circulant_tensor_henry(dx,ko,L,M,N,fl_no_fft);
%             else
%                 fl_no_fft=1;
%                 [fN_all2,st_sparse_precon2] = lse_generate_circulant_tensor_henry(dx,1,L,M,N,fl_no_fft);
% %                 load fN.mat
%                 % note multiply fN_all and st_sparse_precon with ko^2 and compute its FFT
%             end
%         end
%         if num_diels>0
%             ss=zeros(num_diels,1);
%             for ii=1:num_diels
%                 if ii < num_diels
%                     ss(ii)=(dx^2)*((epsa(ii)+epsa(ii+1))/((epsa(ii)-epsa(ii+1))*2*eps0)); % self-term for dielectric matrix
%                 else
%                     ss(ii)=(dx^2)*((epsa(ii)+epsb)/((epsa(ii)-epsb)*2*eps0)); % self-term for dielectric matrix
%                 end
%             end
%         else
%             ss=0;
%         end
%         tend = toc(tini);
%         disp(['Total time for getting circulant tensor ::: ' ,num2str(tend)]);
%         
%     end
%     if fl_Vsour==1
%         if (port_no == num_ports)
%             %         Ae = Ae_original;
%             % Remove rows and columns of Ae corresponding to ground and excitation nodes
%             AL(nodeid_4_grnd,:)=0;
%             AL(nodeid_4_injectcurr,:)=0;
%             Aq(nodeid_4_grnd,:)=0;
%             Aq(nodeid_4_injectcurr,:)=0;
%             %         clear Ae_original
%         end
%     end
% end
% return
% disp(['Done... Precomputing LSE data structures'])
% disp('-----------------------------------------------------')
% % ------------------------------------------------------------------------
% %              Slove LSE
% % -------------------------------------------------------------------------
% disp(['Solving LSEs ...'])
% Y_mat=zeros(num_ports,num_ports,num_freq);
% Y_mat2=zeros(num_ports,num_ports,num_freq);
% Z_mat=zeros(num_ports,num_ports,num_freq);
% R_jL_mat=zeros(num_ports,num_ports,num_freq);
% 
% for freq_no=1:num_freq
%     
%     tinisim = tic;
%     if (num_freq > 1)
%         freq = freq_all(freq_no);
%         %if (freq < 1e6)
%         %    tol = 1e-12;
%         %else
%         %    tol = 1e-8;
%         %end
%     end
%     
%     EMconstants
%     disp('-----------------------------------------------------')
%     disp(['Simulation for frequency : ',num2str(freq),' started! ', 'freq pnt: ',num2str(freq_no), ' / ', num2str(num_freq)])
%     
%     %     rhs_vect(nodeid_4_injectcurr,port_no)=1/(wid/dx*hei/dx);
%     %     rhs_vect(nodeid_4_grnd,port_no)=1/(wid/dx*hei/dx);
%     %     [RHS(N_c+N_out+1:N_c+2*N_out,port_no)] = rhs_vect; % I is source
%     
%     % setting new constitutive parameters for new freq
%     Mr = epsilon_r - 1j*sigma_e/(eo*omega); % permittivity
%     Mc = Mr - 1.0; % susceptibility
%     OneoverMc = 1.0 ./ Mc; % one over susceptibility
%     
%     % circulant tensor for the current frequency
%     if (num_freq > 1)
%         %         ko=1;
%         fN_all = fN_all2*(ko^2);
%         fN_all = fft_operator(fN_all);
%         st_sparse_precon = st_sparse_precon2 * (ko^2);
%         %         fN_all = fN_all2*1i*2*pi*freq*mu;
%         %         fN_all = fft_operator(fN_all);
%         %         st_sparse_precon = st_sparse_precon2*1i*2*pi*freq*mu;
%     end
%     %     for ii=1:3
%     %         for jj=1:3
%     %             fN_ch_all{ii,jj}=fN_ch_all{ii,jj};%/(dx^2);
%     %         end
%     %     end
%     %     fN_ch_all=fN_ch_all*(1/(1i*2*pi*freq));
%     %     st_sparse_precon(4)=fN_ch_all{1}(1,1,1)/(1i*2*pi*freq);%%%%%%%%%%how to extend to diel case
%     for port_no=1:num_ports
%         disp(['Solving for port # ',num2str(port_no), ' ...'])
%         
%         
%         % ------------------------------------------------------------------------
%         %              Assign Excitation and Ground Nodes
%         % ------------------------------------------------------------------------
%         
%         [nodeid_4_grnd,nodeid_4_injectcurr]=lse_assign_exc_grnd_nodes(nodeid_lft,nodeid_rght,nodeid_wlcond,num_ports,port_no);
%         
%         % ------------------------------------------------------------------------
%         %     Solve Linear System of Equations Iteratively
%         % -------------------------------------------------------------------------
%         dim_ind=[N_c Nn];
%         %         dim_ind=[N_c N_in N_out];
%         if (port_no == 1)
%             % prepare the preconditioner
%             tinisim = tic;
%             %             num_diel=0;
%             %             [blk_pre,ids_pre,ord_pre]=implement_preconditioner_ver2nn(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel,freq);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% block diagonal preconditioner %%%%%%
%             %             num_diel=0;
%             %             [blk_pre,ids_pre,ord_pre]=implement_preconditioner_ver2nn(dx,num_vox_in_blk,L,M,N,geom_bndry_panels,Eps_inout,inds_glob,num_diel,freq*dx^4);
%             blk_pre=cell(size(ids_pre,1),1);
%             for ll=1:size(ids_pre,1)
%                 blk_pre{ll}=blk_pre2{ll}*(2*pi*freq*dx^4);
%             end
%             lse_sparse_precon_prepare(dx,freq,OneoverMc,idxS3,st_sparse_precon,nodeid_4_grnd,nodeid_4_injectcurr,AL,Aq,dim_ind,blk_pre,ids_pre,ord_pre,fl_Vsour);
%             disp(['Time for prepare precond ::: ' ,num2str(toc(tinisim))]);
%         end
%         
% % %         Solve the system iteratively
% %                 J=rand(N_c+2*Nn,1);
% %                 fACPU   = lse_matvect_mult(J, fN_all,fN_ch_all, A_ind,dim_ind, OneoverMc, dx, freq, idxS5, nodeid_4_grnd, nodeid_4_injectcurr,1,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,num_diels,ss,blk_pre,ids_pre,ord_pre,fl_Vsour);
% %                 return
%         
%         % Define the handle for matvect
%         fACPU   = @(J)lse_matvect_mult(J, fN_all,fN_ch_all, A_ind,dim_ind, OneoverMc, dx, freq, idxS5, nodeid_4_grnd, nodeid_4_injectcurr,1,inds_glob,diel_pnl_id_locs,diel_pnl_ids,ids_panels,num_diels,ss,blk_pre,ids_pre,ord_pre,fl_Vsour);
%         tini = tic;
%         disp(['Iterative solution started ... '])
%         [rhs_vect_sparse_precon]=lse_sparse_precon_multiply(RHS(:,port_no),AL,Aq,nodeid_4_grnd,nodeid_4_injectcurr,blk_pre,ids_pre,ord_pre);
%         [x, flag, relres, iter, resvec] = pgmres(@(J)fACPU(J), rhs_vect_sparse_precon, inner_it, tol, outer_it);
%         tend = toc(tini);
%         disp(['Total time for iterative solution ::: ' ,num2str(tend)]);
%         disp(['Done... Iterative solution'])
%         
%         if (abs(freq_curr_plot-freq)<1e-12 && port_no == 1)
%             x_backup = x;
%         end
%         % ------------------------------------------------------------------------
%         %     Compute the Currents on Port Nodes and the Column in Ymat
%         % -------------------------------------------------------------------------
%         
%         currs_port_yparams=zeros(num_ports,1);
%         for kk=1:num_ports
%             currs_port_yparams(kk,1)=sum(RHS(:,kk).*x);
%             %             phi_exc=x(N_c+N_out+nodeid_4_injectcurr);
%             %             phi_grnd=x(N_c+N_out+nodeid_4_grnd);
%             %             It_exc=rhs_vect(nodeid_4_injectcurr);
%             %             It_grnd=rhs_vect(nodeid_4_grnd);
%             %             currs_port_yparams(kk,1)=(sum(phi_exc)-sum(phi_grnd))/(wid/dx*hei/dx);%/(It_exc(1));
%         end
%         
%         Y_mat(:,port_no,freq_no)=currs_port_yparams(:,1);
%         
%         % remove the following later on! not needed!
%         % compute column w/ alternative way - just for double checking
%         %         [currs_port_yparams2] = lse_compute_Y_mat_column_alternative(num_ports,Ae,Ae_only_leaving,Ae_only_entering_bndry,x,nodeid_lft,currs_port_yparams);
%         
%         
%         disp(['Done... Solving for port # ',num2str(port_no)])
%         
%     end
%     disp('-----------------------------------------------------')
%     disp(['Done... Simulation for frequency : ',num2str(freq),' freq pnt: ',num2str(freq_no), ' / ', num2str(num_freq)])
%     if fl_Vsour==2 % I source
%         Z_mat(:,:,freq_no)=(squeeze(Y_mat(:,:,freq_no))); % I is source
%     else
%         Z_mat(:,:,freq_no)=inv(squeeze(Y_mat(:,:,freq_no))); %V is source
%     end
%     R_jL_mat(:,:,freq_no)=(real(squeeze(Z_mat(:,:,freq_no))))+sqrt(-1)*(imag(Z_mat(:,:,freq_no)));%/(2*pi*freq));
%     
% end
% RJLMAT=squeeze(R_jL_mat)
% ZZ=sqrt(real(RJLMAT).^2+imag(RJLMAT).^2)
% % ZZ=sqrt(real(R_jL_mat).^2+imag(R_jL_mat).^2);
% % for freq_no=1:num_freq
% %     disp(['R+jL matrix for frequency = ',num2str(freq_all(freq_no))])
% %     for kk=1:num_ports
% %         disp([num2str(R_jL_mat(kk,:,freq_no))])
% %     end
% % end
% 
% 
% 
% 
