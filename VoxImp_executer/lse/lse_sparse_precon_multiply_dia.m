
function [JOut_full_out]=lse_sparse_precon_multiply_dia(JOut_full_in,Ae,Aq)

global A_inv LL UU PP QQ RR Sch_sparse slct_decomp_sch fl_cholmod A_temp

% ---------------------------------------------------------------------
% Sparse preconditioner [E F; G H]
% ---------------------------------------------------------------------

fl_fast_multiply = 1; % for single inversion of Schur complement
fl_volt_source = 2; % symmetric voltage source implementation

num_node=size(Ae,1);
num_curr=size(Ae,2);
JOut_full_out=zeros(num_node+num_curr,1);

switch slct_decomp_sch
    
    case 'lu_decomp'
        
        if (fl_volt_source == 1) % voltage source
            Ae_tmp=Ae; % Attention::: temporarily doubling the memory for Ae!!!
%             Ae_tmp(nodeid_4_grnd,:)=0;
%             Ae_tmp(nodeid_4_injectcurr,:)=0;
            
            if (fl_fast_multiply == 1)
                % for inverting system [A_inv B; C D][x;y]=[a;b]
                % Solve for y of (D - C A_inv B)y = b-C(A_inv)a;
                % Get potentials first. Then obtain currents via
                % (A_inv) x = a-By
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae_tmp * A_inv * JOut_full_in(1:num_curr) ) )))));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    QQ * (UU \ (LL \ (PP * (RR \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    -QQ * (UU \ (LL \ (PP * (RR \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    +QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
            end
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                temp1=Ae * A_inv * JOut_full_in(1:num_curr);
                    temp_unk=JOut_full_in(num_curr+1:num_curr+num_node);
                    for ll=1:size(ids_pre,1)
                        temp_unk(ids_pre{ll},:)=blk_pre{ll}*temp_unk(ids_pre{ll},:);
                    end
                    temp2=Aq*temp_unk;
                    temp=[temp1;temp2];
%                     bb_dum = JOut_full_in(num_curr+num_node+1:num_curr+2*num_node) - temp;
%                     
%                     JOut_full_out(num_curr+num_node+QQ) = ldlsolve (LL,bb_dum(QQ));
                    
                    temp_vect=(JOut_full_in(1:num_curr+num_node) - (([-Ae';-Aq'])*JOut_full_out(num_curr+num_node+1:num_curr+2*num_node)));
                    temp_vect1=temp_vect(1:num_curr).* diag(A_inv);
                    temp_vect2=temp_vect(1+num_curr:num_curr+num_node);
                    for ll=1:size(ids_pre,1)
                        temp_vect2(ids_pre{ll},:)=blk_pre{ll}*temp_vect2(ids_pre{ll},:);
                    end
%                     JOut_full_out(1:num_curr+num_node) = [temp_vect1;temp_vect2];
                    
                    
                    
                JOut_full_out(num_curr+num_node+1:num_curr+2*num_node) = ...
                    QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+num_node+1:num_curr+2*num_node) - temp )))));
                
                JOut_full_out(1:num_curr+num_node) =[temp_vect1;temp_vect2];% (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
                
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    QQ * (UU \ (LL \ (PP * (RR \ (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    -QQ * (UU \ (LL \ (PP * (RR \ (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    +QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
            end
            
        end
        
        
    case 'ldlt_decomp'
        
        if (fl_volt_source == 1) % voltage source
            
            disp('Change decomposition !!! ')
            disp('LDLT decomposition can not work with voltage source !!!')
            error('as the Schur complement is not symmetric for such case ... ')
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                
                if (fl_cholmod == 1)

                    %JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    %    (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) ))))));
                    temp1=Ae * A_inv * JOut_full_in(1:num_curr);
                    temp_unk=JOut_full_in(num_curr+1:num_curr+num_node);
%                     for ll=1:size(ids_pre,1)
%                         temp_unk(ids_pre{ll},:)=blk_pre{ll}*temp_unk(ids_pre{ll},:);
%                     end
                    temp2=A_temp*temp_unk;
                    temp=temp1+temp2;
                    bb_dum = JOut_full_in(num_curr+num_node+1:num_curr+2*num_node) - temp;
                    
                    JOut_full_out(num_curr+num_node+QQ) = ldlsolve (LL,bb_dum(QQ));
                    
                    temp_vect=(JOut_full_in(1:num_curr+num_node) - (([-Ae';-Aq'])*JOut_full_out(num_curr+num_node+1:num_curr+2*num_node)));
                    temp_vect1=temp_vect(1:num_curr).* diag(A_inv);
                    temp_vect2=temp_vect(1+num_curr:num_curr+num_node);
                    temp_vect2=A_temp*temp_vect2;
%                     for ll=1:size(ids_pre,1)
%                         temp_vect2(ids_pre{ll},:)=blk_pre{ll}*temp_vect2(ids_pre{ll},:);
%                     end
                    JOut_full_out(1:num_curr+num_node) = [temp_vect1;temp_vect2];
                else
                    
                    JOut_full_out(num_curr+1:num_curr+num_node) = ...
                        (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) ))))));
                    
                    JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                end
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (Ae*A_inv*JOut_full_in(1:num_curr)))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node)))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * (JOut_full_in(num_curr+1:num_curr+num_node)))))));
                
            end
            
        end
        
    case 'chol_decomp'
        
        if (fl_volt_source == 1)
            
            disp('Change decomposition !!! ')
            disp('Cholesky decomposition can not work with voltage source !!!')
            error('as the Schur complement is not symmetric for such case ... ')
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    (PP * (LL' \ (LL \ (PP' * (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae * A_inv * JOut_full_in(1:num_curr) ) )))));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
            else
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (PP * (LL' \ (LL \ (PP' * (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (PP * (LL' \ (LL \ (PP' * (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - (PP * (LL' \ (LL \ (PP' * (Ae*A_inv*JOut_full_in(1:num_curr))))));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (PP * (LL' \ (LL \ (PP' * (JOut_full_in(num_curr+1:num_curr+num_node))))));
                
            end
        end
        
        
    case 'AGMG'
        
        
        if (fl_volt_source == 1)
            
            Ae_tmp=Ae; % Attention::: temporarily doubling the memory for Ae!!!
%             Ae_tmp(nodeid_4_grnd,:)=0;
%             Ae_tmp(nodeid_4_injectcurr,:)=0;
            
            if (fl_fast_multiply == 1)
                
                disp('Fast solution !!!')
                
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node) - ( Ae_tmp * A_inv * JOut_full_in(1:num_curr) ) ));
                
                JOut_full_out(1:num_curr) = (JOut_full_in(1:num_curr) - ((-Ae')*JOut_full_out(num_curr+1:num_curr+num_node))) .* diag(A_inv);
                
                
            else
                
                disp('Slow solution !!!')
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (Sch_sparse \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr)));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - Sch_sparse \ (Ae_tmp*A_inv*JOut_full_in(1:num_curr));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
            end
            
        else % current source or symmetic voltage source
            
            if (fl_fast_multiply == 1)
                temp1=Ae * A_inv * JOut_full_in(1:num_curr);
                temp_unk=JOut_full_in(num_curr+1:num_curr+num_node);
                temp2=A_temp*temp_unk;
                temp=temp1+temp2;
                bb_dum = JOut_full_in(num_curr+num_node+1:num_curr+2*num_node) - temp;
                [JOut_full_out(num_curr+num_node+1:num_curr+2*num_node),flag,relres,iter,resvec]=agmg((Sch_sparse),bb_dum, ...
                    1,1e-4,100,[],[],202);
                temp_vect=(JOut_full_in(1:num_curr+num_node) - (([-Ae';-Aq'])*JOut_full_out(num_curr+num_node+1:num_curr+2*num_node)));
                temp_vect1=temp_vect(1:num_curr).* diag(A_inv);
                temp_vect2=temp_vect(1+num_curr:num_curr+num_node);
                temp_vect2=A_temp*temp_vect2;
                JOut_full_out(1:num_curr+num_node) = [temp_vect1;temp_vect2];
                
            else
                
                %disp('Slow solution !!!')
                
                % block E contribution
                JOut_full_out(1:num_curr) = A_inv*JOut_full_in(1:num_curr)+A_inv*(-Ae')*...
                    (Sch_sparse \ (Ae*A_inv*JOut_full_in(1:num_curr)));
                
                % block F contribution
                JOut_full_out(1:num_curr) = JOut_full_out(1:num_curr)...
                    + A_inv * (Ae') * (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
                % block G contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = ...
                    - Sch_sparse \ (Ae*A_inv*JOut_full_in(1:num_curr));
                
                % block H contribution
                JOut_full_out(num_curr+1:num_curr+num_node) = JOut_full_out(num_curr+1:num_curr+num_node)...
                    + (Sch_sparse \ (JOut_full_in(num_curr+1:num_curr+num_node)));
                
            end
        end
    case 'no_decomp'
        temp1=Ae * A_inv * JOut_full_in(1:num_curr);
        temp_unk=JOut_full_in(num_curr+1:num_curr+num_node);
        for ll=1:size(ids_pre,1)
            temp_unk(ids_pre{ll},:)=blk_pre{ll}*temp_unk(ids_pre{ll},:);
        end
        %                     temp2=Aq*temp_unk;
        temp2=temp_unk;
        temp=temp1+temp2;
        bb_dum = JOut_full_in(num_curr+num_node+1:num_curr+2*num_node) - temp;
        
        a_temp=strumpack_solve(real(bb_dum));
        b_temp=strumpack_solve(imag(bb_dum));
        JOut_full_out(num_curr+num_node+1:num_curr+2*num_node)=a_temp+sqrt(-1).*b_temp;
        
        temp_vect=(JOut_full_in(1:num_curr+num_node) - (([-Ae';-Aq'])*JOut_full_out(num_curr+num_node+1:num_curr+2*num_node)));
        temp_vect1=temp_vect(1:num_curr).* diag(A_inv);
        temp_vect2=temp_vect(1+num_curr:num_curr+num_node);
        for ll=1:size(ids_pre,1)
            temp_vect2(ids_pre{ll},:)=blk_pre{ll}*temp_vect2(ids_pre{ll},:);
        end
        JOut_full_out(1:num_curr+num_node) = [temp_vect1;temp_vect2];
    otherwise
        
        error('Invalid decomposition selection for Schur complement')
        
end

fprintf ('.') ;
