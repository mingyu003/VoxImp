function [bndry_panel]=port_geom_panels(dx,grid_intcon,idxS,pnt_rght,pnt_lft,port_direction)

% the function is used to generate geometry panels and boundary panels
% bndry_panel: col 1-3: coordinates, col 4: normal direction, col 5:
% numbering, col 6: # of conductor, col 7: # of dielectric

[L,M,N,~,~]=size(grid_intcon);
% num_conds=max(idxS(:,2));
boolean_tens=zeros(L,M,N);
boolean_tens(idxS(:,1))=1;
% %%%---------------------------------------------------------------
% %                    find port id
% %%%---------------------------------------------------------------
pnt_rght=cell2mat(pnt_rght);
pnt_lft=cell2mat(pnt_lft);
size_lft=size(pnt_lft,1);
size_rght=size(pnt_rght,1);

bndry_panel=zeros(size_lft+size_rght,10);

for ii=1: size_lft
    bndry_panel(ii,1:3)=pnt_lft(ii,:);
    bndry_panel(ii,4)=port_direction(1);
    bndry_panel(ii,6)=1; % need to be generalized
end
for ii=1:size_rght
    bndry_panel(size_lft+ii,1:3)=pnt_rght(ii,:);
    bndry_panel(size_lft+ii,4)=port_direction(2);
    bndry_panel(size_lft+ii,6)=1; % need to be generalized
end
    
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
 
    
    
    
    
    
    
% disp('Done... Extracting unique panels of a grid')
% disp('-----------------------------------------------------')