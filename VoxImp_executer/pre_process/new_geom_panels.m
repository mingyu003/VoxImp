function [bndry_panel]=new_geom_panels(dx,grid_intcon,idxS)

% the function is used to generate geometry panels and boundary panels
% bndry_panel: col 1-3: coordinates, col 4: normal direction, col 5:
% numbering, col 6: # of conductor, col 7: # of dielectric

[L,M,N,~,~]=size(grid_intcon);
% num_conds=max(idxS(:,2));
%%%------------------------------------------------------------------
%             obtain geom panels' coordinates and normal
%%%------------------------------------------------------------------
boolean_tens=zeros(L,M,N);
boolean_tens(idxS(:,1))=1;

dumxr=0; % x right
dumyb=0; % y back
dumzu=0; % z up
dumxl=0; % x left
dumyf=0; % y front
dumzb=0; % z bottom

bndry_x_r=zeros(L*M*N,10);
bndry_x_l=zeros(L*M*N,10);
bndry_y_f=zeros(L*M*N,10);
bndry_y_b=zeros(L*M*N,10);
bndry_z_b=zeros(L*M*N,10);
bndry_z_u=zeros(L*M*N,10);
total_panels= total_panel(boolean_tens, grid_intcon,dx);
%total_panels has 30 cols for the 6 faces of a voxel, every 5 cols for one
%face. 1-3 cols mean the coordinates, 4 and 5 cols mean the cond and diel
%numbering

for kk=1:L
    for ll=1:M
        for mm=1:N
            if (boolean_tens(kk,ll,mm))==1%-1)<1e-12
                
                if ((kk+1)<=L && boolean_tens(kk+1,ll,mm)<1e-12) || kk==L  % right bndry
                    dumxr=dumxr+1;
                    bndry_x_r(dumxr,1:3)=total_panels(kk,ll,mm,6:8);
                    bndry_x_r(dumxr,4)=1;
                    bndry_x_r(dumxr,6)=total_panels(kk,ll,mm,9);
                    bndry_x_r(dumxr,7)=total_panels(kk,ll,mm,10);
                end

                if ((ll+1)<=M &&  boolean_tens(kk,ll+1,mm)<1e-12) || ll==M  %back bndry  
                    dumyb=dumyb+1;
                    bndry_y_b(dumyb,1:3)=total_panels(kk,ll,mm,16:18);
                    bndry_y_b(dumyb,4)=2;
                    bndry_y_b(dumyb,6)=total_panels(kk,ll,mm,19);
                    bndry_y_b(dumyb,7)=total_panels(kk,ll,mm,20);
                end
                
                if ((mm+1)<=N && boolean_tens(kk,ll,mm+1)<1e-12) || mm==N%up bndry
                    dumzu=dumzu+1;
                    bndry_z_u(dumzu,1:3)=total_panels(kk,ll,mm,26:28);
                    bndry_z_u(dumzu,4)=3;
                    bndry_z_u(dumzu,6)=total_panels(kk,ll,mm,29);
                    bndry_z_u(dumzu,7)=total_panels(kk,ll,mm,30);
                end
                
                if ((kk-1)>=1 && boolean_tens(kk-1,ll,mm)<1e-12) || kk==1 %left bndry
                    dumxl=dumxl+1;
                    bndry_x_l(dumxl,1:3)=total_panels(kk,ll,mm,1:3);
                    bndry_x_l(dumxl,4)=1;
                    bndry_x_l(dumxl,6)=total_panels(kk,ll,mm,4);
                    bndry_x_l(dumxl,7)=total_panels(kk,ll,mm,5);
                end
                
                if ((ll-1)>=1 && boolean_tens(kk,ll-1,mm)<1e-12)|| ll==1%front bndry
                    dumyf=dumyf+1;
                    bndry_y_f(dumyf,1:3)=total_panels(kk,ll,mm,11:13);
                    bndry_y_f(dumyf,4)=2;
                    bndry_y_f(dumyf,6)=total_panels(kk,ll,mm,14);
                    bndry_y_f(dumyf,7)=total_panels(kk,ll,mm,15);
                end
                if ((mm-1)>=1 && boolean_tens(kk,ll,mm-1)<1e-12) || mm==1 %bottom bndry
                    dumzb=dumzb+1;
                    bndry_z_b(dumzb,1:3)=total_panels(kk,ll,mm,21:23);
                    bndry_z_b(dumzb,4)=3;
                    bndry_z_b(dumzb,6)=total_panels(kk,ll,mm,24);
                    bndry_z_b(dumzb,7)=total_panels(kk,ll,mm,25);
                end
            end
        end  
    end
end

bndry_x_r=bndry_x_r(1:dumxr,:);
bndry_x_l=bndry_x_l(1:dumxl,:);
bndry_y_f=bndry_y_f(1:dumyf,:);
bndry_y_b=bndry_y_b(1:dumyb,:);
bndry_z_b=bndry_z_b(1:dumzb,:);
bndry_z_u=bndry_z_u(1:dumzu,:);

bndry_panel=[bndry_x_r; bndry_x_l; bndry_y_b; bndry_y_f; bndry_z_u; bndry_z_b];

%%%-------------------------------------------------------------------
%                  obtain numbering for geom panels
%%%-------------------------------------------------------------------
[tensor_ijk_xdir_pnl,tensor_ijk_ydir_pnl,tensor_ijk_zdir_pnl]=panel_numbering(L,M,N);
bbox_origin=[0,0,0];
org_x_panel=[bbox_origin(1) bbox_origin(2)+dx/2 bbox_origin(3)+dx/2];
org_y_panel=[bbox_origin(1)+dx/2 bbox_origin(2) bbox_origin(3)+dx/2];
org_z_panel=[bbox_origin(1)+dx/2 bbox_origin(2)+dx/2 bbox_origin(3)];

%%%---------------------------------------------------------------
%                    find bundry panels numbering in computational domain
%%%---------------------------------------------------------------
for kk=1:size(bndry_panel,1)
    if (bndry_panel(kk,4) == 1) % x-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_x_panel)/dx)+1;
        bndry_panel(kk,5)=tensor_ijk_xdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (bndry_panel(kk,4) == 2) % y-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_y_panel)/dx)+1;
        bndry_panel(kk,5)=tensor_ijk_ydir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (bndry_panel(kk,4) == 3) % z-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_z_panel)/dx)+1;
        bndry_panel(kk,5)=tensor_ijk_zdir_pnl(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    end
end
%%%---------------------------------------------------------------
%                    find bundry panels numbering in geom domain
%%%---------------------------------------------------------------
[tensor_ijk_xdir_pnl_g,tensor_ijk_ydir_pnl_g,tensor_ijk_zdir_pnl_g]=panel_numbering_geom(L,M,N,boolean_tens);

for kk=1:size(bndry_panel,1)
    if (bndry_panel(kk,4) == 1) % x-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_x_panel)/dx)+1;
        bndry_panel(kk,10)=tensor_ijk_xdir_pnl_g(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (bndry_panel(kk,4) == 2) % y-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_y_panel)/dx)+1;
        bndry_panel(kk,10)=tensor_ijk_ydir_pnl_g(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    elseif (bndry_panel(kk,4) == 3) % z-directed panel
        tmp_ind=round((bndry_panel(kk,1:3)-org_z_panel)/dx)+1;
        bndry_panel(kk,10)=tensor_ijk_zdir_pnl_g(tmp_ind(1),tmp_ind(2),tmp_ind(3));
    end
end
% %%%---------------------------------------------------------------
% %                    find port id
% %%%---------------------------------------------------------------
% no_port=size(pnt_rght,1);
% for ii=no_port
%     size_rght=size(pnt_rght{ii},1);
%     size_lft=size(pnt_lft{ii},1);
%     for jj=1:size_rght
%         
%     
%     
%     
    
    
    
    
    
    
    
    
% disp('Done... Extracting unique panels of a grid')
% disp('-----------------------------------------------------')