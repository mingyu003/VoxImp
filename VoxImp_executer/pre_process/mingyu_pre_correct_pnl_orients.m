function [geom_bndry_panels] =mingyu_pre_correct_pnl_orients(num_cond,num_diel,geom_bndry_panels,Cnt,dx,fl_plot_pnl_orients)
% for ii=1:num_cond
%     %     geom_bndry_panels{ii}(:,7)=0;
%     geom_bndry_panels{ii}(:,8)=Eps_inout(1,1);
%     geom_bndry_panels{ii}(:,9)=Eps_inout(1,2);
% end
% if num_diel>0
%     for ii=1:num_diel
%         geom_bndry_panels{num_cond+ii}(:,7)=ii;
%         geom_bndry_panels{num_cond+ii}(:,8)=Eps_inout(num_cond+ii,1);
%         geom_bndry_panels{num_cond+ii}(:,9)=Eps_inout(num_cond+ii,2);
%     end
% end
for ii=1:num_diel
    ref_cen=Cnt(num_cond+ii,1:3);
    for ll=1:size(geom_bndry_panels{num_cond+ii},1)
        cen_pnl=geom_bndry_panels{num_cond+ii}(ll,1:3);
        orient_pnl=geom_bndry_panels{num_cond+ii}(ll,4);
        vect_ref2pnl=(cen_pnl-ref_cen)/norm(cen_pnl-ref_cen);
        
        if (orient_pnl == 1) % x-directed panel
            unit_vect=[1 0 0];
        elseif (orient_pnl == 2) % y-directed panel
            unit_vect=[0 1 0];
        elseif (orient_pnl == 3) % z-directed panel
            unit_vect=[0 0 1];
        end
        
        ang_btw_vects=acos(sum(vect_ref2pnl.*unit_vect)/(norm(vect_ref2pnl)*norm(unit_vect)));
        
        if (ang_btw_vects > pi/2) % change the sign
            geom_bndry_panels{num_cond+ii}(ll,4)=-geom_bndry_panels{num_cond+ii}(ll,4);
        end
    end
end




if (fl_plot_pnl_orients == 1)
    
    geom_pnl_info=cell2mat(geom_bndry_panels);
    for ll=1:size(geom_pnl_info,1)
        if (ll == 1)
            figure
            set(gca,'FontSize',24)
        end
        
        if (orient_pnl == 1) % x-directed panel
            unit_vect=[1 0 0];
        elseif (orient_pnl == 2) % y-directed panel
            unit_vect=[0 1 0];
        elseif (orient_pnl == 3) % z-directed panel
            unit_vect=[0 0 1];
        end
        
        if (sign(geom_pnl_info(ll,4))==-1)
            unit_vect=-unit_vect;
        end
        cen_pnl=geom_pnl_info(ll,1:3);
        %if (cen_pnl(1,1) > 0.5) % uncomment this if you'd like to cut fig
        quiver3(cen_pnl(1),cen_pnl(2),cen_pnl(3),unit_vect(1)*dx,unit_vect(2)*dx,unit_vect(3)*dx)
        hold on
        %end
        
        % lets put the panel as well
        
        if (orient_pnl == 1) % x-directed
            verts=[cen_pnl(1) cen_pnl(2)-dx/2 cen_pnl(3)-dx/2; ...
                cen_pnl(1) cen_pnl(2)+dx/2 cen_pnl(3)-dx/2; ...
                cen_pnl(1) cen_pnl(2)+dx/2 cen_pnl(3)+dx/2; ...
                cen_pnl(1) cen_pnl(2)-dx/2 cen_pnl(3)+dx/2;];
        elseif (orient_pnl == 2) % y-directed
            verts=[cen_pnl(1)-dx/2 cen_pnl(2) cen_pnl(3)-dx/2; ...
                cen_pnl(1)+dx/2 cen_pnl(2) cen_pnl(3)-dx/2; ...
                cen_pnl(1)+dx/2 cen_pnl(2) cen_pnl(3)+dx/2; ...
                cen_pnl(1)-dx/2 cen_pnl(2) cen_pnl(3)+dx/2;];
        elseif (orient_pnl == 3) % z-directed
            verts=[cen_pnl(1)-dx/2 cen_pnl(2)-dx/2 cen_pnl(3); ...
                cen_pnl(1)+dx/2 cen_pnl(2)-dx/2 cen_pnl(3); ...
                cen_pnl(1)+dx/2 cen_pnl(2)+dx/2 cen_pnl(3); ...
                cen_pnl(1)-dx/2 cen_pnl(2)+dx/2 cen_pnl(3);];
        end
        
        %if (verts(1,1) > 0.5 && verts(2,1) > 0.5 && verts(3,1) > 0.5 && verts(4,1) > 0.5) % uncomment this if you'd like to cut fig
        patch(verts(:,1),verts(:,2),verts(:,3),[1 2 4 3],'EdgeColor','blue','FaceColor','none');
        hold on
        %end
        
        if (ll == size(geom_pnl_info,1))
            xlabel('x');ylabel('y');zlabel('z');
            title('Surface Normals of Panels');
            grid on
            set(gca,'FontSize',24)
            axis tight
        end
    end
    
end

