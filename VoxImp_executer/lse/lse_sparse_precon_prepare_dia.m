function lse_sparse_precon_prepare_dia(dx,freq,OneoverMc,idxS3,st_sparse_precon,Ae,dim_ind,dia_cap)
global A_inv LL UU PP QQ RR Sch_sparse slct_decomp_sch fl_cholmod A_temp

%slct_decomp_sch='ldlt_decomp'; %'chol_decomp','no_decomp','ldlt_decomp','chol_decomp'

% slct_decomp_sch='no_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp'
slct_decomp_sch='ldlt_decomp'; %'no_decomp','lu_decomp','ldlt_decomp','chol_decomp','AGMG'
fl_cholmod = 1; % flag for using CHOLMOD in suitesparse for fast and memory-efficient LDLT factorization and inversion
fl_volt_source = 2; % symmetric voltage source implementation
fl_profile = 0; % cpu and memory profiling


fl_accu_check = 0; % flag for accuracy checking for inversion with the decompositions
tic_Assembly = tic;

% constants
% N_c = dim_ind(1); %current unknowns % 5*non_air_voxel
% N_in = dim_ind(2); % # internal nodes
% N_out = dim_ind(3); % # external nodes
Nc = dim_ind(1); %current unknowns % 5*non_air_voxel
Nn = dim_ind(2); % # internal nodes

omega = 2*pi*freq;
mu = 4*pi*1e-7;
co = 299792458;
eo = 1/co^2/mu;

% 1) Compute A_inv matrix
tic
OneoverMc_dum=OneoverMc(idxS3(1));
diag_pulse=1/((1/(j*omega*eo))*((1/dx*OneoverMc_dum)-st_sparse_precon(1)));
diag_2Dlinear=1/((1/(j*omega*eo))*((dx/6*OneoverMc_dum*(1/(dx^2)))-st_sparse_precon(2)));
diag_3Dlinear=1/((1/(j*omega*eo))*((dx/2*OneoverMc_dum*(1/(dx^2)))-st_sparse_precon(3)));

inds=zeros(Nc,3);
inds(1:Nc,1)=[1:1:Nc];
inds(1:Nc,2)=inds(1:Nc,1);
inds(1:3*Nc/5,3)=abs(diag_pulse);
inds(3*Nc/5+1:4*Nc/5,3)=abs(diag_2Dlinear);
inds(4*Nc/5+1:Nc,3)=abs(diag_3Dlinear);
A_inv=sparse(inds(:,1),inds(:,2),inds(:,3));

if (fl_profile == 1); disp(['Time to compute A_inv block ::: ',num2str(toc)]); end

infomem1 = whos('A_inv');
memestimated = (infomem1.bytes)/(1024*1024);
if (fl_profile == 1); disp(['Memory for A_inv block (MB)::' , num2str(memestimated)]); end

%%%%%%%%%%%%% blk diagonal
% nz=0;
% for ll=1:size(ids_pre,1)
%     nz=nz+length(ids_pre{ll})*length(ids_pre{ll});
% end
% A_temp=spalloc(Nn,Nn,nz);
ind1=1:Nn; ind2=1:Nn; 
aa=dia_cap*ones(Nn,1);
A_temp=sparse(ind1,ind2,aa);
% transpose_Aq=Aq';
% for ll=1:size(ids_pre,1)
%     transpose_Aq(ids_pre{ll},ids_pre{ll})=blk_pre{ll}*transpose_Aq(ids_pre{ll},ids_pre{ll});
% end

%     if isempty(ord_pre)==0
%         for kk=1:size(ord_pre,1)
%             ch_coef_dum(ord{kk},:)=(1/ss(end))*ch_coef_dum(ord{kk},:);
%         end
%     end

% 2) Obtain Schur complement
tic
% if (fl_volt_source == 1)
%     DD = spalloc(Nn,Nn,length(nodeid_4_grnd)+length(nodeid_4_injectcurr));
% else
%     DD = spalloc(Nn,Nn,length(nodeid_4_grnd));
% end

% for kk=1:length(nodeid_4_grnd)
%     DD(nodeid_4_grnd(kk),nodeid_4_grnd(kk))=1;
% end

% if (fl_volt_source == 1 || fl_volt_source == 2)
%     for kk=1:length(nodeid_4_injectcurr)
%         DD(nodeid_4_injectcurr(kk),nodeid_4_injectcurr(kk))=1;
%     end
% end

% if (fl_volt_source == 1)
%     A_ind_tmp=Ae; % Attention::: temporarily doubling the memory for A_ind!!!
%     A_ind_tmp(nodeid_4_grnd,:)=0;
%     A_ind_tmp(nodeid_4_injectcurr,:)=0;
%     Sch_comp=DD+A_ind_tmp*A_inv*Ae';
%     clear A_ind_tmp
% else % symm. voltage source and current source will use this
    Sch_comp1=Ae*A_inv*Ae';
    Sch_comp2=A_temp;
    Sch_comp=Sch_comp1+Sch_comp2;
% end
if (fl_profile == 1); disp(['Time to obtain Schur complement ::: ',num2str(toc)]); end


% 3) Factorization of Schur complement with different schemes

%ishermitian(Sch_comp)
if (issymmetric(Sch_comp) == 0)
    if (strcmp(slct_decomp_sch, 'ldlt_decomp') == 1)
        disp('Schur complement matrix is not symmetric')
        disp('Decomposition is changed to LU')
        slct_decomp_sch = 'lu_decomp';
    elseif (strcmp(slct_decomp_sch, 'chol_decomp') == 1 )
        disp('Schur complement matrix is not symmetric')
        disp('Decomposition is changed to LU')
        slct_decomp_sch = 'lu_decomp';
    end
end

switch slct_decomp_sch
    
    case 'lu_decomp'
        
        tic
        [LL,UU,PP,QQ,RR] = lu(Sch_comp);
        disp(['Time to (LU) factorize Schur complement ::: ',num2str(toc)])
        
        infomem1 = whos('LL');infomem2 = whos('UU'); infomem3 = whos('PP');
        infomem4 = whos('QQ');infomem5 = whos('RR');
        memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes+infomem4.bytes+infomem5.bytes)/(1024*1024);
        disp(['Memory for LU fact. matrices (MB)::' , num2str(memestimated)]);
        % Test for CPU time of inversion
        bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
        tic
        x_dum_lu = QQ * (UU \ (LL \ (PP * (RR \ bb))));
        disp(['Time for one inversion with LU ::: ',num2str(toc)])
        tic
        
        if (fl_accu_check == 1)
            tic
            x_dum_back = Sch_comp\bb;
            disp(['Time for one inversion with backslash ::: ',num2str(toc)])
            csd=num2str(max(abs(x_dum_lu-x_dum_back)./abs(x_dum_back)));
            disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
        end
        
        
    case 'ldlt_decomp'
        
        if (fl_volt_source == 1)
            error('Sch comp not symmetric for voltage source - LDLT decomp does not work - change decomposition')
        end
        
        if (fl_cholmod == 1)
            % factorization and inversion with CHOLMOD in Suitesparse
            tic;
            [LL,p,QQ] = ldlchol (Sch_comp);
            disp(['Time to (LDLT) factorize Schur complement - Cholmod ::: ',num2str(toc)])
            
            infomem1 = whos('LL');infomem2 = whos('p');infomem3 = whos('QQ');
            memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes)/(1024*1024);
            disp(['Memory for LDLT fact. matrices - Cholmod (MB)::' , num2str(memestimated)]);
            
            %x_dum_ldlt=zeros(size(Sch_comp,1),1);
            %bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
            
            %tic;
            %x_dum_ldlt(QQ) = ldlsolve (LL,bb(QQ));
            %disp(['Time for one inversion with LDLT - Cholmod ::: ',num2str(toc)])
            
            if (fl_accu_check == 1)
                tic
                x_dum_back = Sch_comp\bb;
                disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                csd=num2str(max(abs(x_dum_ldlt-x_dum_back)./abs(x_dum_back)));
                disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
            end
            
        else
            
            tic
            [LL,QQ,PP,RR] = ldl(Sch_comp);% option 2 PP'*RR*A*RR*PP = LL*QQ*LL'
            disp(['Time to (LDLT) factorize Schur complement ::: ',num2str(toc)])
            
            infomem1 = whos('LL');infomem2 = whos('QQ');infomem3 = whos('PP'); infomem4 = whos('RR');
            memestimated = (infomem1.bytes+infomem2.bytes+infomem3.bytes+infomem4.bytes)/(1024*1024);
            disp(['Memory for LDLT fact. matrices (MB)::' , num2str(memestimated)]);
            % Test for CPU time of inversion
            bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
            tic
            x_dum_ldlt = RR * PP * (LL' \ (QQ \ (LL \ (PP' * RR * bb)))); % for option 2
            disp(['Time for one inversion with LDLT ::: ',num2str(toc)])
            
            if (fl_accu_check == 1)
                tic
                x_dum_back = Sch_comp\bb;
                disp(['Time for one inversion with backslash ::: ',num2str(toc)])
                csd=num2str(max(abs(x_dum_ldlt-x_dum_back)./abs(x_dum_back)));
                disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
            end
        end
        
        
    case 'chol_decomp'
        
        if (fl_volt_source == 1)
            error('Sch comp not symmetric for voltage source - Cholesky decomp does not work - change decomposition')
        end
        
        tic
        [LL,g,PP] = chol(Sch_comp,'lower'); % option 3: LL'=S'AS
        if (g ~= 0)
            error ('MATLAB:posdef', 'Matrix must be positive definite.') ;
        end
        disp(['Time to (Cholesky) factorize Schur complement ::: ',num2str(toc)])
        
        infomem1 = whos('LL');infomem2 = whos('PP');
        memestimated = (infomem1.bytes+infomem2.bytes)/(1024*1024);
        disp(['Memory for Chol. fact. matrices (MB)::' , num2str(memestimated)]);
        % Test for CPU time of inversion
        %bb=randn(size(Sch_comp,1),1)+sqrt(-1)*randn(size(Sch_comp,1),1);
        %tic
        %x_dum_chol = PP * (LL' \ (LL \ (PP' * bb))); % for option 3
        %disp(['Time for one inversion with Chol ::: ',num2str(toc)])
        
        if (fl_accu_check == 1)
            tic
            x_dum_back = Sch_comp\bb;
            disp(['Time for one inversion with backslash ::: ',num2str(toc)])
            csd=num2str(max(abs(x_dum_chol-x_dum_back)./abs(x_dum_back)));
            disp(['Max rel diff between solutions obtained via decomposition and direct baskslash ::: ',csd])
        end
        
        % For large scA_ind test, the vector permutation implementation below
        % didn't reduce the time. However, I include this here for future
        % reference
        
        % % Cholesky with vector permutation
        % tic
        % [LL,g,s] = chol(Sch_comp,'lower','vector');
        % disp(['Time to (Cholesky) factorize Schur complement - vector ::: ',num2str(toc)])
        %
        % x_dum_chol3=zeros(size(x_dum_chol,1),size(x_dum_chol,2));
        %
        % tic
        % x_dum_chol3(s) = (LL' \ (LL \ (bb(s)))); % for option 3
        % disp(['Time for one inversion with Chol - vector ::: ',num2str(toc)])
        % tic
        %
        % rel_diff2=max(abs(x_dum_chol -x_dum_chol3)./abs(x_dum_chol))
        
        
    case 'AGMG'
        
        infomem1 = whos('Sch_comp');
        memestimated = (infomem1.bytes)/(1024*1024);
        disp(['Memory for Schur complement matrix (MB)::' , num2str(memestimated)]);
        Sch_sparse = Sch_comp;
                Sch_sparse=spfun(@(x)complex(x,1e-8),Sch_sparse);
                agmg((Sch_sparse),[],1,[],[],[],[],1);
     case 'no_decomp'
%         spparms('spumoni',2)
        infomem1 = whos('Sch_comp');
        memestimated = (infomem1.bytes)/(1024*1024);
        disp(['Memory for Schur complement matrix (MB)::' , num2str(memestimated)]);
        b_temp=rand(Nn,1);%+sqrt(-1)*rand(Nn,1);
        x_temp=strumpack_solve(Sch_comp,b_temp);
    otherwise
        
        error('Invalid decomposition selection for Schur complement')
        
end

Time_Assembly = toc(tic_Assembly);
disp(['Time to prepare sparse precon ::: ',num2str(Time_Assembly)])