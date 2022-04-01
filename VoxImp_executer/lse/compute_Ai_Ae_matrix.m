function [Ae, Ai,nodeid_lft,nodeid_rght,nodeid_wlcond,Ai_only_leaving,Ai_only_entering_bndry] = compute_Ai_Ae_matrix(idxS, grid_intcon, L, M, N, dx, pnt_lft, pnt_rght,port_direct,index_bndry_node,fl_check_ports)

% % Additional codelet for seperate running
% clc; close all; clear all; 
% L=3; M=2; N=2; dx = 1;
% % L=8; M=3; N=4; dx = 1;
% grid_tmp = ones(L,M,N);
% % grid_tmp(4:5,:,3:4)=0; 
% idxS = find(abs(grid_tmp(:)) > 1e-12); clear grid_tmp; 
% port_direct=[1 1];
% [grid_intcon] = generategridfrombbox(dx,[0 L*dx],[0 M*dx],[0 N*dx],0);
% fl_meth_comp=2 % methods: 0-> slow , 1-> fast with sort, 2-> fast w/ graph
% pnt_lft={[0 0.5 0.5; 0 1.5 0.5; 0 0.5 1.5; 0 1.5 1.5]};pnt_well_cond=[];
% pnt_rght={[3 0.5 0.5; 3 1.5 0.5; 3 0.5 1.5; 3 1.5 1.5]};fl_check_ports=0;


nodeid_4_injectcurr=[];
nodeid_4_grnd=[];
nodeid_4_well_cond=[];

tstart = tic;
fl_profile = 0;

% constants
num_nonair_cube=length(idxS);

% Inputs for visualization
fl_vis_only_voxels_with_ids=0;
fl_vis_panel_locs_ids_x=0;
fl_vis_panel_locs_ids_y=0;
fl_vis_panel_locs_ids_z=0;
fl_vis_exc_grnd = fl_check_ports;

num_elem_allowed = 5000;
freq_vis = 1;
if (num_nonair_cube > num_elem_allowed)% visualization of boxes is slow when nD > num_elem_allowed
    freq_vis = 25; % every xx element unknown id will be written.
end

% putting port nodes in temporary arrays

% find number of ports
num_ports=size(pnt_lft,1);
if (size(pnt_rght,1) ~= size(pnt_lft,1))
    error('Port definition is wrong !!! Number of faces should 2 for each port')
end

% find number of nodes on each face of each port
num_node_each_port=zeros(num_ports,2);
for kk=1:num_ports
    num_node_each_port(kk,1)=size(pnt_rght{kk},1);
    num_node_each_port(kk,2)=size(pnt_lft{kk},1);
end

% temporarily define the right nodes as ground
pnt_grnd=zeros(sum(num_node_each_port(:,1)),3); % right points
% temporarily define the left nodes as excitation
pnt_exc=zeros(sum(num_node_each_port(:,2)),3); % left point

% put the nodes in temporary arrays
dum_exc=0;
dum_grnd=0;
for kk=1:num_ports
    pnt_exc(dum_exc+1:dum_exc+num_node_each_port(kk,2),1:3) = pnt_lft{kk}(:,1:3);
    dum_exc=dum_exc+num_node_each_port(kk,2);
    
    pnt_grnd(dum_grnd+1:dum_grnd+num_node_each_port(kk,1),1:3) = pnt_rght{kk}(:,1:3);
    dum_grnd=dum_grnd+num_node_each_port(kk,1);
end

%     % 1) Get the boolean tensor
boolean_tens=zeros(L,M,N);
boolean_tens(idxS)=1;

dumxr_all=0; % x right
dumyb_all=0; % y back
dumzu_all=0; % z up
dumxl_all=0; % x left
dumyf_all=0; % y front
dumzb_all=0; % z bottom

all_x_r=zeros(L*M*N,5);
all_x_l=zeros(L*M*N,5);
all_y_f=zeros(L*M*N,5);
all_y_b=zeros(L*M*N,5);
all_z_b=zeros(L*M*N,5);
all_z_u=zeros(L*M*N,5);
total_panels= all_panel(boolean_tens, grid_intcon,dx);
dumxr=0; dumyb=0; dumzu=0; dumxl=0; dumyf=0; dumzb=0;
non_air_vox_ind=0;
for mm=1:N
    for ll=1:M
        for kk=1:L
            if (boolean_tens(kk,ll,mm))>1e-12%-1)<1e-12
                non_air_vox_ind=non_air_vox_ind+1;
                
                dumxr_all=dumxr_all+1;
                all_x_r(dumxr_all,1:3)=total_panels(kk,ll,mm,5:7);
                all_x_r(dumxr_all,4)=1;
                if ((kk+1)<=L && boolean_tens(kk+1,ll,mm)<1e-12) || kk==L  % right bndry
                    dumxr=dumxr+1;
                end
                
                dumyb_all=dumyb_all+1;
                all_y_b(dumyb_all,1:3)=total_panels(kk,ll,mm,13:15);
                all_y_b(dumyb_all,4)=2;
                if ((ll+1)<=M &&  boolean_tens(kk,ll+1,mm)<1e-12) || ll==M  %back bndry
                    dumyb=dumyb+1;
                end
                dumzu_all=dumzu_all+1;
                all_z_u(dumzu_all,1:3)=total_panels(kk,ll,mm,21:23);
                all_z_u(dumzu_all,4)=3;
                if ((mm+1)<=N && boolean_tens(kk,ll,mm+1)<1e-12) || mm==N%up bndry
                    dumzu=dumzu+1;
                end
                dumxl_all=dumxl_all+1;
                all_x_l(dumxl_all,1:3)=total_panels(kk,ll,mm,1:3);
                all_x_l(dumxl_all,4)=1;
                if ((kk-1)>=1 && boolean_tens(kk-1,ll,mm)<1e-12) || kk==1 %left bndry
                    dumxl=dumxl+1;
                end
                dumyf_all=dumyf_all+1;
                all_y_f(dumyf_all,1:3)=total_panels(kk,ll,mm,9:11);
                all_y_f(dumyf_all,4)=2;
                if ((ll-1)>=1 && boolean_tens(kk,ll-1,mm)<1e-12)|| ll==1%front bndry
                    dumyf=dumyf+1;
                end
                dumzb_all=dumzb_all+1;
                all_z_b(dumzb_all,1:3)=total_panels(kk,ll,mm,17:19);
                all_z_b(dumzb_all,4)=3;
                if ((mm-1)>=1 && boolean_tens(kk,ll,mm-1)<1e-12) || mm==1 %bottom bndry
                    dumzb=dumzb+1;
                end
            end
        end
    end
end

all_x_r=all_x_r(1:dumxr_all,:);
all_x_l=all_x_l(1:dumxl_all,:);
all_y_f=all_y_f(1:dumyf_all,:);
all_y_b=all_y_b(1:dumyb_all,:);
all_z_b=all_z_b(1:dumzb_all,:);
all_z_u=all_z_u(1:dumzu_all,:);

all_panels=[all_x_r; all_x_l; all_y_b; all_y_f; all_z_u; all_z_b];

%%%-------------------------------------------------------------------
%                  obtain numbering for geom panels
%%%-------------------------------------------------------------------
tensor_ijk_xdir_pnl=zeros(L+1,M,N);
dum_cnt=1;
for mm=1:N% z-variation
    for ll=1:M % y-variation
        for kk=1:L+1 % x-variation
            if kk<=L
                if (boolean_tens(kk,ll,mm))>1e-12
                    tensor_ijk_xdir_pnl(kk,ll,mm)=dum_cnt;
                    dum_cnt=dum_cnt+1;
                else
                    if kk>1
                        if (boolean_tens(kk-1,ll,mm))>1e-12
                            tensor_ijk_xdir_pnl(kk,ll,mm)=dum_cnt;
                            dum_cnt=dum_cnt+1;
                        end
                    end
                end
            elseif kk==L+1 && (boolean_tens(kk-1,ll,mm))>1e-12
                tensor_ijk_xdir_pnl(kk,ll,mm)=dum_cnt;
                dum_cnt=dum_cnt+1;
            end
        end
    end
end
% generate ids of  y_directed panels
tensor_ijk_ydir_pnl=zeros(L,M+1,N);
for mm=1:N% z-variation
    for ll=1:M+1 % y-variation
        for kk=1:L % x-variation
% for kk=1:L % x-variation 
%     for mm=1:N % z-variation 
%         for ll=1:M+1 % y-variation
            if ll<=M
                if (boolean_tens(kk,ll,mm))>1e-12
                    tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
                    dum_cnt=dum_cnt+1;
                else
                    if ll>1
                        if (boolean_tens(kk,ll-1,mm))>1e-12
                            tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
                            dum_cnt=dum_cnt+1;
                        end
                    end
                end
            elseif ll==M+1 && (boolean_tens(kk,ll-1,mm))>1e-12
                tensor_ijk_ydir_pnl(kk,ll,mm)=dum_cnt;
                dum_cnt=dum_cnt+1;
            end
        end
    end
end

% generate ids of z_directed panels
tensor_ijk_zdir_pnl=zeros(L,M,N+1);
for mm=1:N+1 % z-variation
    for ll=1:M % y-variation
        for kk=1:L % x-variation
% for ll=1:M % y-variation 
%     for kk=1:L % x-variation 
%         for mm=1:N+1 % z-variation
            if mm<=N
                if (boolean_tens(kk,ll,mm))>1e-12 && mm<=N
                    tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
                    dum_cnt=dum_cnt+1;
                else
                    if mm>1
                        if (boolean_tens(kk,ll,mm-1))>1e-12
                            tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
                            dum_cnt=dum_cnt+1;
                        end
                    end
                end
            elseif mm==N+1 && (boolean_tens(kk,ll,mm-1))>1e-12
                tensor_ijk_zdir_pnl(kk,ll,mm)=dum_cnt;
                dum_cnt=dum_cnt+1;
            end
        end
    end
end

bbox_origin=[0,0,0];
org_x_panel=[bbox_origin(1) bbox_origin(2)+dx/2 bbox_origin(3)+dx/2];
org_y_panel=[bbox_origin(1)+dx/2 bbox_origin(2) bbox_origin(3)+dx/2];
org_z_panel=[bbox_origin(1)+dx/2 bbox_origin(2)+dx/2 bbox_origin(3)];

for kk=1:size(all_panels,1)
    if (all_panels(kk,4) == 1) % x-directed panel
        tmp_ind=round((all_panels(kk,1:3)-org_x_panel)/dx)+1;
        all_panels(kk,5)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (all_panels(kk,4) == 2) % y-directed panel
        tmp_ind=round((all_panels(kk,1:3)-org_y_panel)/dx)+1;
        all_panels(kk,5)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (all_panels(kk,4) == 3) % z-directed panel
        tmp_ind=round((all_panels(kk,1:3)-org_z_panel)/dx)+1;
        all_panels(kk,5)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%           port nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nodeid_4_injectcurr=zeros(size(pnt_exc,1),1);
for kk=1:size(pnt_exc,1)
    if (port_direct(1) == 1) % x-directed panel
        tmp_ind=round((pnt_exc(kk,1:3)-org_x_panel)/dx)+1;
        nodeid_4_injectcurr(kk)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (port_direct(1) == 2) % y-directed panel
        tmp_ind=round((pnt_exc(kk,1:3)-org_y_panel)/dx)+1;
        nodeid_4_injectcurr(kk)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (port_direct(1) == 3) % z-directed panel
        tmp_ind=round((pnt_exc(kk,1:3)-org_z_panel)/dx)+1;
        nodeid_4_injectcurr(kk)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    end
end

nodeid_4_grnd=zeros(size(pnt_grnd,1),1);
for kk=1:size(pnt_grnd,1)
    if (port_direct(2) == 1) % x-directed panel
        tmp_ind=round((pnt_grnd(kk,1:3)-org_x_panel)/dx)+1;
        nodeid_4_grnd(kk)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (port_direct(2) == 2) % y-directed panel
        tmp_ind=round((pnt_grnd(kk,1:3)-org_y_panel)/dx)+1;
        nodeid_4_grnd(kk)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (port_direct(2) == 3) % z-directed panel
        tmp_ind=round((pnt_grnd(kk,1:3)-org_z_panel)/dx)+1;
        nodeid_4_grnd(kk)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voxel2panel_x=zeros(dumxr_all,2);
voxel2panel_y=zeros(dumyb_all,2);
voxel2panel_z=zeros(dumzu_all,2);
voxel2panel_x(:,2)=all_panels(1:dumxr_all,5);
voxel2panel_x(:,1)=all_panels(dumxr_all+1:dumxr_all+dumxl_all,5);
voxel2panel_y(:,2)=all_panels(dumxr_all+dumxl_all+1:dumxr_all+dumxl_all+dumyb_all,5);
voxel2panel_y(:,1)=all_panels(dumxr_all+dumxl_all+dumyb_all+1:dumxr_all+dumxl_all+dumyb_all+dumyf_all,5);
voxel2panel_z(:,2)=all_panels(dumxr_all+dumxl_all+dumyb_all+dumyf_all+1:dumxr_all+dumxl_all+dumyb_all+dumyf_all+dumzu_all,5);
voxel2panel_z(:,1)=all_panels(dumxr_all+dumxl_all+dumyb_all+dumyf_all+dumzu_all+1:dumxr_all+dumxl_all+dumyb_all+dumyf_all+dumzu_all+dumzb_all,5);

all_panels_ids=[voxel2panel_x; voxel2panel_y; voxel2panel_z];
if(fl_profile == 1); disp(['Time for indexing panels ::: ', num2str(toc)]); end;

% 8) forming Ai matrix
% num_panel_x=max(max(voxel2panel_x));
% num_panel_y=max(max(voxel2panel_y));
% num_panel_z=max(max(voxel2panel_z));
tic
num_nodes=max(max(voxel2panel_z));%num_panel_x+num_panel_y+num_panel_z;

% enter +1, leaves -1; first entry for leaving, second entry for entering
const_lin=1/2;
%const_lin=dx/2;

sp_mat_inds=zeros(16*num_nonair_cube,3);
sp_mat_inds_only_leaving_currs=zeros(8*num_nonair_cube,3);
for kk=1:3*num_nonair_cube % pertinent to Jx, Jy, and Jz
    sp_mat_inds(kk,1:3)= [all_panels_ids(kk,1) kk -1];
    sp_mat_inds(3*num_nonair_cube+kk,1:3)= [all_panels_ids(kk,2) kk 1];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(kk,1:3)= [all_panels_ids(kk,1) kk -1];
end

dum=1;
for kk=3*num_nonair_cube+1:4*num_nonair_cube % pertinent to Jfx
    sp_mat_inds(6*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
    sp_mat_inds(7*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(3*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
    dum=dum+1;
end

dum=1;
for kk=3*num_nonair_cube+1:4*num_nonair_cube % pertinent to Jfy
    sp_mat_inds(8*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
    sp_mat_inds(9*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk -const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(4*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk -const_lin];
    dum=dum+1;
end

% for space diagonal currents
dum=1;
for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsx
    sp_mat_inds(10*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
    sp_mat_inds(11*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,2) kk const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(5*num_nonair_cube+dum,1:3)= [all_panels_ids(dum,1) kk const_lin];
    dum=dum+1;
end

dum=1;
for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsy
    sp_mat_inds(12*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk const_lin];
    sp_mat_inds(13*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,2) kk const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(6*num_nonair_cube+dum,1:3)= [all_panels_ids(num_nonair_cube+dum,1) kk const_lin];
    dum=dum+1;
end

dum=1;
for kk=4*num_nonair_cube+1:5*num_nonair_cube % pertinent to Jsz
    sp_mat_inds(14*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,1) kk -2*const_lin];
    sp_mat_inds(15*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,2) kk -2*const_lin];
    % the following is added for current vis.
    sp_mat_inds_only_leaving_currs(7*num_nonair_cube+dum,1:3)= [all_panels_ids(2*num_nonair_cube+dum,1) kk -2*const_lin];
    dum=dum+1;
end

%Ai=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3));
Ai=sparse(sp_mat_inds(:,1),sp_mat_inds(:,2),sp_mat_inds(:,3),num_nodes,5*num_nonair_cube);

num_bndry_nodes=dumxr+dumyb+dumzu+dumxl+dumyf+dumzb;
% Ae=sparse(num_bndry_nodes,5*num_nonair_cube);
% counter=0;jj=0;
Ae=Ai(index_bndry_node(:,1),:);
Ai(index_bndry_node(:,1),:)=[];
jj=1;
for ii=1:num_nodes
    if ii==index_bndry_node(jj,1) %boundary nodes
%         counter=counter+1;
        if isempty(find(nodeid_4_injectcurr==ii))==0
            nodeid_4_injectcurr(find(nodeid_4_injectcurr==ii))=jj;
        end
        if isempty(find(nodeid_4_grnd==ii))==0
            nodeid_4_grnd(find(nodeid_4_grnd==ii))=jj;
        end
        if jj<length(index_bndry_node)
            jj=jj+1;
        end
    end
end
% for ii=1:num_nodes
%     jj=jj+1;
%     if abs(sum(Ai(jj,1:3*num_nonair_cube)))>1e-10 %boundary nodes
%         counter=counter+1;
%         if isempty(find(nodeid_4_injectcurr==ii))==0
%             nodeid_4_injectcurr(find(nodeid_4_injectcurr==ii))=counter;
%         end
%         if isempty(find(nodeid_4_grnd==ii))==0
%             nodeid_4_grnd(find(nodeid_4_grnd==ii))=counter;
%         end
%         Ae(counter,:)=Ai(jj,:);
%         Ai(jj,:)=[];
%         jj=jj-1;
%     end
% end
%in Ai, col 1:3*num_nonair_voxel, if the non-zero element is in k-th
%column, it means k is the number of this boundary nodes. Note that the
%order of Phi_ne should be same as this.

infomem1 = whos('Ai');
memestimated = (infomem1.bytes)/(1024*1024);
disp(['Memory for Ai matrix (MB)::' , num2str(memestimated)]);

% Additional data structure for visualization of currents

% We need two data structures:
% (1) Ai matrix for leaving currents
Ai_only_leaving=sparse(sp_mat_inds_only_leaving_currs(:,1),sp_mat_inds_only_leaving_currs(:,2),sp_mat_inds_only_leaving_currs(:,3),num_nodes,5*num_nonair_cube);

% (2) Ai matrix for entering currents to boundary nodes
% Finding the nodes to which current enters
dum_vect=zeros(size(Ai,2),1);
dum_vect(1:size(Ai,2)/5*3)=1;
dum_res=Ai*dum_vect; % ids of boundary nodes -1:exiting, +1:entering

nodes_only_curr_enter=find(dum_res == 1);
sp_mat_dum=Ai(nodes_only_curr_enter,:);
[rows,cols,vals] = find(sp_mat_dum);
for kk=1:length(rows) %correct rows
    rows(kk)=nodes_only_curr_enter(rows(kk));
end

Ai_only_entering_bndry=sparse(rows,cols,vals,num_nodes,5*num_nonair_cube);


% The following is for visualization!!!
% retrieving the centers of non-air voxels
tic
dum=1;
xy_curr_ids_locs=zeros(3*num_nonair_cube,3);
for mm=1:N
    for ll=1:M
        for kk=1:L
            if (grid_intcon(kk,ll,mm,1) ~=0 && grid_intcon(kk,ll,mm,2) ~=0 && ...
                    grid_intcon(kk,ll,mm,3) ~=0 ) % This is dangerous!
                coor_tmp=squeeze(grid_intcon(kk,ll,mm,1:3));
                xy_curr_ids_locs(dum,1:3)=coor_tmp;
                xy_curr_ids_locs(num_nonair_cube+dum,1:3)=coor_tmp;
                xy_curr_ids_locs(2*num_nonair_cube+dum,1:3)=coor_tmp;
                dum=dum+1;
            end
        end
    end
end

% finding locations of panels enclosing x and y currents
tic
all_panel_locs=zeros(3*num_nonair_cube,6);
for kk=1:num_nonair_cube
    %Jx currents
    all_panel_locs(kk,1:3)=[xy_curr_ids_locs(kk,1)-dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
    all_panel_locs(kk,4:6)=[xy_curr_ids_locs(kk,1)+dx*0.5 xy_curr_ids_locs(kk,2) xy_curr_ids_locs(kk,3)];
    %Jy currents
    ind=num_nonair_cube+kk;
    all_panel_locs(ind,1:3)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)-dx*0.5 xy_curr_ids_locs(ind,3)];
    all_panel_locs(ind,4:6)=[xy_curr_ids_locs(ind,1) xy_curr_ids_locs(ind,2)+dx*0.5 xy_curr_ids_locs(ind,3)];
    %Jz currents
    ind2=2*num_nonair_cube+kk;
    all_panel_locs(ind2,1:3)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)-dx*0.5];
    all_panel_locs(ind2,4:6)=[xy_curr_ids_locs(ind2,1) xy_curr_ids_locs(ind2,2) xy_curr_ids_locs(ind2,3)+dx*0.5];
end
if(fl_profile == 1); disp(['Time for finding locations of surface panels ::: ', num2str(toc)]); end

% visualize graphs along x, y, or z directions
fl_vis_graphs=0;
if(fl_vis_graphs == 1)
    figure; plot(G_x);title('Graph for panels along x direction')
    figure; plot(G_y);title('Graph for panels along y direction')
    figure; plot(G_z);title('Graph for panels along z direction')
end

if(fl_profile == 1); disp(['Total time for generating Ai matrix ::: ', num2str(toc(tstart))]); end;


%% Finding the nodes on which excitation and ground points are defined

disp('-----------------------------------------------------')
disp('Finding the IDs of port nodes... ')


% putting the node ids of ports nodes in structure

% puting node/current ids of each port in a structure
nodeid_lft=cell(num_ports,1);
nodeid_rght=cell(num_ports,1);

dum_exc=0;
dum_grnd=0;
for kk=1:num_ports
    nodeid_lft{kk}(:) = nodeid_4_injectcurr(dum_exc+1:dum_exc+num_node_each_port(kk,2));
    dum_exc=dum_exc+num_node_each_port(kk,2);
    
    nodeid_rght{kk}(:) = nodeid_4_grnd(dum_grnd+1:dum_grnd+num_node_each_port(kk,1));
    dum_grnd=dum_grnd+num_node_each_port(kk,1);
end

nodeid_wlcond = nodeid_4_well_cond;

disp('Done... Finding the IDs of port nodes... ')
disp('-----------------------------------------------------')

%if (fl_check_ports == 1); error('Just plotting port nodes... No Simulation...'); end

num_node = size(Ai,1);
num_curr = size(Ai,2);
disp('-----------------------------------------------------')
disp(['Number of current unknowns ::: ',num2str(num_curr)])
disp(['Number of potential unknowns ::: ',num2str(num_node)])
disp(['Number of total unknowns ::: ',num2str(num_curr+num_node)])
disp('-----------------------------------------------------')
