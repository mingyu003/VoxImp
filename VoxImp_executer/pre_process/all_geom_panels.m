function [all_panels]=all_geom_panels(dx,grid_intcon,idxS,num_vox_in_blk)

% the function is used to generate geometry panels and boundary panels
% bndry_panel: col 1-3: coordinates, col 4: normal direction, col 5:
% numbering, col 6: # of conductor, col 7: # of dielectric

[L,M,N,~,~]=size(grid_intcon);
num_blk_x=ceil(L/num_vox_in_blk);  num_blk_y=ceil(M/num_vox_in_blk);  num_blk_z=ceil(N/num_vox_in_blk);
bound=floor([L,M,N]./num_vox_in_blk)*num_vox_in_blk;

if mod(L,num_vox_in_blk)>0 && L>1
    num_big_blk_x=num_blk_x-1;
else
    num_big_blk_x=num_blk_x;
end
if mod(M,num_vox_in_blk)>0 && M>1
    num_big_blk_y=num_blk_y-1;
else
    num_big_blk_y=num_blk_y;
end
if mod(N,num_vox_in_blk)>0 && N>1
    num_big_blk_z=num_blk_z-1;
else
    num_big_blk_z=num_blk_z;
end

num_blk_x=num_big_blk_x;  num_blk_y=num_big_blk_y;  num_blk_z=num_big_blk_z;
% % num_conds=max(idxS(:,2));
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

bndry_x_r=zeros(L*M*N,11);
bndry_x_l=zeros(L*M*N,11);
bndry_y_f=zeros(L*M*N,11);
bndry_y_b=zeros(L*M*N,11);
bndry_z_b=zeros(L*M*N,11);
bndry_z_u=zeros(L*M*N,11);
%col 11 indicates the panel is boundary or not
total_panels= total_panel(boolean_tens, grid_intcon,dx);
%total_panels has 30 cols for the 6 faces of a voxel, every 5 cols for one
%face. 1-3 cols mean the coordinates, 4 and 5 cols mean the cond and diel
%numbering

for kk=1:L
    for ll=1:M
        for mm=1:N
            if (boolean_tens(kk,ll,mm))==1%-1)<1e-12
                dumxl=dumxl+1;
                bndry_x_l(dumxl,1:3)=total_panels(kk,ll,mm,1:3);
                bndry_x_l(dumxl,4)=1;
                bndry_x_l(dumxl,6)=total_panels(kk,ll,mm,4);
                bndry_x_l(dumxl,7)=total_panels(kk,ll,mm,5);
                if kk>bound(1)
                    temp(1)=floor(kk/num_vox_in_blk);
                else
                    temp(1)=ceil(kk/num_vox_in_blk);
                end
                if ll>bound(2)
                    temp(2)=floor(ll/num_vox_in_blk);
                else
                    temp(2)=ceil(ll/num_vox_in_blk);
                end
                if mm>bound(3)
                    temp(3)=floor(mm/num_vox_in_blk);
                else
                    temp(3)=ceil(mm/num_vox_in_blk);
                end
                %                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                bndry_x_l(dumxl,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                if ((kk-1)>=1 && boolean_tens(kk-1,ll,mm)<1e-12) || kk==1 %left bndry
                    bndry_x_l(dumxl,11)=1; % boundary panels
                end
                if ((kk+1)<=L && boolean_tens(kk+1,ll,mm)<1e-12) || kk==L  % right bndry
                    dumxr=dumxr+1;
                    bndry_x_r(dumxr,1:3)=total_panels(kk,ll,mm,6:8);
                    bndry_x_r(dumxr,4)=1;
                    bndry_x_r(dumxr,6)=total_panels(kk,ll,mm,9);
                    bndry_x_r(dumxr,7)=total_panels(kk,ll,mm,10);
                    bndry_x_r(dumxr,11)=1;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
                    %                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_x_r(dumxr,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                
                dumyf=dumyf+1;
                bndry_y_f(dumyf,1:3)=total_panels(kk,ll,mm,11:13);
                bndry_y_f(dumyf,4)=2;
                bndry_y_f(dumyf,6)=total_panels(kk,ll,mm,14);
                bndry_y_f(dumyf,7)=total_panels(kk,ll,mm,15);
                if kk>bound(1)
                    temp(1)=floor(kk/num_vox_in_blk);
                else
                    temp(1)=ceil(kk/num_vox_in_blk);
                end
                if ll>bound(2)
                    temp(2)=floor(ll/num_vox_in_blk);
                else
                    temp(2)=ceil(ll/num_vox_in_blk);
                end
                if mm>bound(3)
                    temp(3)=floor(mm/num_vox_in_blk);
                else
                    temp(3)=ceil(mm/num_vox_in_blk);
                end
                %                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                bndry_y_f(dumyf,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                if ((ll-1)>=1 && boolean_tens(kk,ll-1,mm)<1e-12)|| ll==1%front bndry
                    bndry_y_f(dumyf,11)=1;
                end
                if ((ll+1)<=M &&  boolean_tens(kk,ll+1,mm)<1e-12) || ll==M  %back bndry
                    dumyb=dumyb+1;
                    bndry_y_b(dumyb,1:3)=total_panels(kk,ll,mm,16:18);
                    bndry_y_b(dumyb,4)=2;
                    bndry_y_b(dumyb,6)=total_panels(kk,ll,mm,19);
                    bndry_y_b(dumyb,7)=total_panels(kk,ll,mm,20);
                    bndry_y_b(dumyb,11)=1;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
                    %                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_y_b(dumyb,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                end
                
                dumzb=dumzb+1;
                bndry_z_b(dumzb,1:3)=total_panels(kk,ll,mm,21:23);
                bndry_z_b(dumzb,4)=3;
                bndry_z_b(dumzb,6)=total_panels(kk,ll,mm,24);
                bndry_z_b(dumzb,7)=total_panels(kk,ll,mm,25);
                if kk>bound(1)
                    temp(1)=floor(kk/num_vox_in_blk);
                else
                    temp(1)=ceil(kk/num_vox_in_blk);
                end
                if ll>bound(2)
                    temp(2)=floor(ll/num_vox_in_blk);
                else
                    temp(2)=ceil(ll/num_vox_in_blk);
                end
                if mm>bound(3)
                    temp(3)=floor(mm/num_vox_in_blk);
                else
                    temp(3)=ceil(mm/num_vox_in_blk);
                end
                %                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                bndry_z_b(dumzb,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
                if ((mm-1)>=1 && boolean_tens(kk,ll,mm-1)<1e-12) || mm==1 %bottom bndry
                    bndry_z_b(dumzb,11)=1;
                end
                if ((mm+1)<=N && boolean_tens(kk,ll,mm+1)<1e-12) || mm==N%up bndry
                    dumzu=dumzu+1;
                    bndry_z_u(dumzu,1:3)=total_panels(kk,ll,mm,26:28);
                    bndry_z_u(dumzu,4)=3;
                    bndry_z_u(dumzu,6)=total_panels(kk,ll,mm,29);
                    bndry_z_u(dumzu,7)=total_panels(kk,ll,mm,30);
                    bndry_z_u(dumzu,11)=1;
                    if kk>bound(1)
                        temp(1)=floor(kk/num_vox_in_blk);
                    else
                        temp(1)=ceil(kk/num_vox_in_blk);
                    end
                    if ll>bound(2)
                        temp(2)=floor(ll/num_vox_in_blk);
                    else
                        temp(2)=ceil(ll/num_vox_in_blk);
                    end
                    if mm>bound(3)
                        temp(3)=floor(mm/num_vox_in_blk);
                    else
                        temp(3)=ceil(mm/num_vox_in_blk);
                    end
                    %                    temp=ceil([kk,ll,mm]/num_vox_in_blk);
                    bndry_z_u(dumzu,10)=temp(1)+(temp(2)-1)*num_blk_x+(temp(3)-1)*num_blk_x*num_blk_y;
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

all_panels=[bndry_x_r; bndry_x_l; bndry_y_b; bndry_y_f; bndry_z_u; bndry_z_b];

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
% %%%---------------------------------------------------------------
% %                    find bundry panels numbering in geom domain
% %%%---------------------------------------------------------------
% [tensor_ijk_xdir_pnl_g,tensor_ijk_ydir_pnl_g,tensor_ijk_zdir_pnl_g]=panel_numbering_geom(L,M,N,boolean_tens);
% 
% for kk=1:size(all_panels,1)
%     if (all_panels(kk,4) == 1) % x-directed panel
%         tmp_ind=round((all_panels(kk,1:3)-org_x_panel)/dx)+1;
%         all_panels(kk,10)=tensor_ijk_xdir_pnl_g(tmp_ind(1),tmp_ind(2),tmp_ind(3));
%     elseif (all_panels(kk,4) == 2) % y-directed panel
%         tmp_ind=round((all_panels(kk,1:3)-org_y_panel)/dx)+1;
%         all_panels(kk,10)=tensor_ijk_ydir_pnl_g(tmp_ind(1),tmp_ind(2),tmp_ind(3));
%     elseif (all_panels(kk,4) == 3) % z-directed panel
%         tmp_ind=round((all_panels(kk,1:3)-org_z_panel)/dx)+1;
%         all_panels(kk,10)=tensor_ijk_zdir_pnl_g(tmp_ind(1),tmp_ind(2),tmp_ind(3));
%     end
% end
% 

% disp('Done... Extracting unique panels of a grid')
% disp('-----------------------------------------------------')