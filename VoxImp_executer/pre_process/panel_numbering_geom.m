function [tensor_ijk_xdir_pnl,tensor_ijk_ydir_pnl,tensor_ijk_zdir_pnl]=panel_numbering_geom(L,M,N,boolean_tens)

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
