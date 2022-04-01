function [JOut_full]=lse_matvect_mult_ver0(JIn0, fN_all,fN_ch_all, A_ind,dim_ind, OneoverMc, dx, freq, idx, nodeid_4_grnd,nodeid_4_injectcurr,fl_sparseprecon,inds_glob,locs_diel_pos,ids_panel_diel,ids_panels,num_diel,ss,blk_pre,ids_pre,ord_pre)


% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------

fl_volt_source = 2; %symmetric voltage source
% fl_sparseprecon = 0;
fl_gpu = 0;
fl_profile = 0;

tic
% constants
mu = 4*pi*1e-7; co = 299792458; eo = 1/co^2/mu; omega = 2*pi*freq;
oneoverjomegaeo=1/(1j*omega*eo);
jomegamu=1j*omega*mu;
oneoversigma_e=OneoverMc*oneoverjomegaeo;

N_c=dim_ind(1); N_in=dim_ind(2);% N_out=dim_ind(3);
% JIn0(N_c+N_out+1:N_c+2*N_out)=JIn0(N_c+N_out+1:N_c+2*N_out)*dx^2;
% fft dimensions
[LfN, MfN, NfN, ~] = size(fN_all{1});
for ii=1:7
    fN_all{ii}=fN_all{ii}*jomegamu;
end
% [LfN, MfN, NfN, ~] = size(fN_all);
[LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
[LfN_xy, MfN_xy, NfN_xy, ~] = size(fN_ch_all{1,2,1}); %xy
[LfN_xz, MfN_xz, NfN_xz, ~] = size(fN_ch_all{1,3,1}); %xz
% domain dimensions
[L, M, N] = size(OneoverMc);

JIn = zeros(L, M, N, 5);
JOut = zeros(L, M, N, 5);
JIn_x = zeros(L+1, M, N);

if (num_diel > 0)
    ch_coefs_xx = JIn0(N_c+inds_glob(1,1):N_c+inds_glob(2,2)); % [rhoa_x_c;rhoa_x_d]
    ch_coefs_yy = JIn0(N_c+inds_glob(3,1):N_c+inds_glob(4,2)); % [rhoa_y_c;rhoa_y_d]
    ch_coefs_zz = JIn0(N_c+inds_glob(5,1):N_c+inds_glob(6,2)); % [rhoa_z_c;rhoa_z_d]
    coe_xx=ss*JIn0(N_c+inds_glob(2,1):N_c+inds_glob(2,2)); % coefficient for self-term extraction
    coe_yy=ss*JIn0(N_c+inds_glob(4,1):N_c+inds_glob(4,2));
    coe_zz=ss*JIn0(N_c+inds_glob(6,1):N_c+inds_glob(6,2));
else
    ch_coefs_xx = JIn0(N_c+inds_glob(1,1):N_c+inds_glob(1,2)); % [rhoa_x_c]
    ch_coefs_yy = JIn0(N_c+inds_glob(2,1):N_c+inds_glob(2,2)); % [rhoa_y_c]
    ch_coefs_zz = JIn0(N_c+inds_glob(3,1):N_c+inds_glob(3,2)); % [rhoa_z_c]
end
JIn_x([ids_panels{1}{1};ids_panels{2}{1}]) = ch_coefs_xx;
fJ_x = fftn(JIn_x(:,:,:),[LfN_x, MfN_x, NfN_x]); % fJ_yx = fJ_zx =fJ_x;
clear JIn_x;
JIn_y = zeros(L, M+1, N);
JIn_y([ids_panels{1}{2};ids_panels{2}{2}]) = ch_coefs_yy;
fJ_y = fftn(JIn_y(:,:,:),[LfN_xy, MfN_xy, NfN_xy]); % fJ_zy = fJ_xy = fJ_y;
clear JIn_y;
JIn_z = zeros(L, M, N+1);
JIn_z([ids_panels{1}{3};ids_panels{2}{3}]) = ch_coefs_zz;
fJ_z = fftn(JIn_z(:,:,:),[LfN_xz, MfN_xz, NfN_xz]); % fJ_yz = fJ_xz = fJ_z;
clear JIn_z;

% translate from local (idx) to global (L,M,N) coordinates
JIn(idx) = JIn0(1:N_c);


% JOut_full = zeros(N_c+N_in+2*N_out,1);
% JOut_full_in = zeros(N_c+N_in+2*N_out,1);
JOut_full = zeros(N_c+2*N_in,1);
JOut_full_in = zeros(N_c+2*N_in,1);


%%% allocate space
%JIn = zeros(L, M, N, 3);
%JOut = zeros(L, M, N, 3);
%JOut_full = zeros(N_c+N_in,1);

% translate from local (idx) to global (L,M,N) coordinates
%JIn(idx) = JIn0(1:N_c);


% fN_all=fN_all*jomegamu;
% ---------------------------------------------------------------------
% apply fft and mv-op for each of the components of JIn
% ---------------------------------------------------------------------

% % x component of JIn, store contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
% Jout1 = fN_all(:,:,:,1) .* fJ; % Gxx*Jx
% Jout4 = fN_all(:,:,:,4) .* fJ; % G2dx*Jx
% %Jout5 = fN_all(:,:,:,4) .* fJ; % G3dx*Jx
% Jout5 = Jout4; % G3dx*Jx
% 
% % y component of JIn, add contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
% Jout2 = fN_all(:,:,:,1) .* fJ; % Gyy*Jy
% Jout4 = Jout4 + fN_all(:,:,:,5) .* fJ; % G2dy*Jy - (-(y-yc))
% Jout5 = Jout5 - fN_all(:,:,:,5) .* fJ; % G3dy*Jy - (y-yc)
% 
% % z component of JIn, store contribution on 2 components of Jout
% fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
% Jout3 = fN_all(:,:,:,1) .* fJ; % Gzz*Jz
% Jout5 = Jout5 + fN_all(:,:,:,7) .* fJ; % G3dz*Jz
% 
% % 2d component of JIn, add contribution on 4 components of Jout
% fJ = fftn(JIn(:,:,:,4),[LfN, MfN, NfN]);
% Jout1 = Jout1 + fN_all(:,:,:,2) .* fJ; % Gx2d*J2d
% Jout2 = Jout2 + fN_all(:,:,:,3) .* fJ; % Gy2d*J2d - (-(y-yc))
% Jout4 = Jout4 + fN_all(:,:,:,6) .* fJ; % G2d2d*J2d
% Jout5 = Jout5 + fN_all(:,:,:,9) .* fJ; % G3d2d*J2d
% 
% % 3d component of JIn, add contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,5),[LfN, MfN, NfN]);
% Jout1 = Jout1 + fN_all(:,:,:,2) .* fJ; % Gx3d*J3d
% Jout2 = Jout2 - fN_all(:,:,:,3) .* fJ; % Gy3d*J3d - (y-yc)
% Jout3 = Jout3 + fN_all(:,:,:,8) .* fJ; % Gz3d*J3d
% Jout4 = Jout4 + fN_all(:,:,:,9) .* fJ; % G2d3d*J3d
% Jout5 = Jout5 + fN_all(:,:,:,10) .* fJ; % G3d3d*J3d
%%%%%%%%%%%%%%%%%% Enrico code
% x component of JIn, store contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
Jout1 = fN_all{1} .* fJ; % Gxx*Jx
Jout4 = -fN_all{2} .* fJ; % G2dx*Jx
Jout5 = Jout4; % G3dx*Jx

% y component of JIn, add contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
Jout2 = fN_all{1} .* fJ; % Gyy*Jy
Jout4 = Jout4 - fN_all{3} .* fJ; % G2dy*Jy - (-(y-yc))
Jout5 = Jout5 + fN_all{3} .* fJ; % G3dy*Jy - (y-yc)

% z component of JIn, store contribution on 2 components of Jout
fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
Jout3 = fN_all{1} .* fJ; % Gzz*Jz
Jout5 = Jout5 - fN_all{5} .* fJ; % G3dz*Jz

% 2d component of JIn, add contribution on 4 components of Jout
fJ = fftn(JIn(:,:,:,4),[LfN, MfN, NfN]);
Jout1 = Jout1 + fN_all{2} .* fJ; % Gx2d*J2d
Jout2 = Jout2 + fN_all{3} .* fJ; % Gy2d*J2d - (-(y-yc))
Jout4 = Jout4 + fN_all{4} .* fJ; % G2d2d*J2d
Jout5 = Jout5 + fN_all{6} .* fJ; % G3d2d*J2d

% 3d component of JIn, add contribution on 3 components of Jout
fJ = fftn(JIn(:,:,:,5),[LfN, MfN, NfN]);
Jout1 = Jout1 + fN_all{2} .* fJ; % Gx3d*J3d
Jout2 = Jout2 - fN_all{3} .* fJ; % Gy3d*J3d - (y-yc)
Jout3 = Jout3 + fN_all{5} .* fJ; % Gz3d*J3d
Jout4 = Jout4 + fN_all{6} .* fJ; % G2d3d*J3d
Jout5 = Jout5 + fN_all{7} .* fJ; % G3d3d*J3d

% apply ifft 
Jout1 = ifftn(Jout1);
Jout2 = ifftn(Jout2);
Jout3 = ifftn(Jout3);
Jout4 = ifftn(Jout4);
Jout5 = ifftn(Jout5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
JOut(:,:,:,1) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,1) + Jout1(1:L,1:M,1:N);
JOut(:,:,:,2) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,2) + Jout2(1:L,1:M,1:N);
JOut(:,:,:,3) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,3) + Jout3(1:L,1:M,1:N);
JOut(:,:,:,4) = (1/(6*dx)) .* oneoversigma_e .* JIn(:,:,:,4) + Jout4(1:L,1:M,1:N);
JOut(:,:,:,5) = (1/(2*dx)) .* oneoversigma_e .* JIn(:,:,:,5) + Jout5(1:L,1:M,1:N);

% % apply ifft and add identity term
% Jout1 = ifftn(Jout1);
% JOut(:,:,:,1) = (1/dx) .* OneoverMc .* JIn(:,:,:,1) - Jout1(1:L,1:M,1:N);
% Jout2 = ifftn(Jout2);
% JOut(:,:,:,2) = (1/dx) .* OneoverMc .* JIn(:,:,:,2) - Jout2(1:L,1:M,1:N);
% Jout3 = ifftn(Jout3);
% JOut(:,:,:,3) = (1/dx) .* OneoverMc .* JIn(:,:,:,3) - Jout3(1:L,1:M,1:N);
% Jout4 = ifftn(Jout4);
% JOut(:,:,:,4) = (dx/6/(dx^2)) .* OneoverMc .* JIn(:,:,:,4) - Jout4(1:L,1:M,1:N);
% Jout5 = ifftn(Jout5);
% JOut(:,:,:,5) = (dx/2/(dx^2)) .* OneoverMc .* JIn(:,:,:,5) - Jout5(1:L,1:M,1:N);
% % multiply by 1/(jweps0)
% 
% JOut(:,:,:,:) = JOut(:,:,:,:) * oneoverjomegaeo;
% 
% % Jout1 = ifftn(Jout1);
% % JOut(:,:,:,1) = (1/dx) * oneoversigma_e(1,1,1) .* JIn(:,:,:,1) + Jout1(1:L,1:M,1:N);
% % Jout2 = ifftn(Jout2);
% % JOut(:,:,:,2) = (1/dx) * oneoversigma_e(1,1,1) .* JIn(:,:,:,2) + Jout2(1:L,1:M,1:N);
% % Jout3 = ifftn(Jout3);
% % JOut(:,:,:,3) = (1/dx) * oneoversigma_e(1,1,1) .* JIn(:,:,:,3) + Jout3(1:L,1:M,1:N);
% % Jout4 = ifftn(Jout4);
% % JOut(:,:,:,4) = (dx/6/(dx^2)) * oneoversigma_e(1,1,1) .* JIn(:,:,:,4) + Jout4(1:L,1:M,1:N);
% % Jout5 = ifftn(Jout5);
% % JOut(:,:,:,5) = (dx/2/(dx^2)) * oneoversigma_e(1,1,1) .* JIn(:,:,:,5) + Jout5(1:L,1:M,1:N);


Jout1_x = ifftn(fJ_x .* fN_ch_all{1,1,1}+fJ_y .* fN_ch_all{1,2,1}+fJ_z .* fN_ch_all{1,3,1});% xx,xy,xz
Jout1_x = Jout1_x(1:L+1,1:M,1:N);
JOut_full(N_c+inds_glob(1,1):N_c+inds_glob(1,2)) = (Jout1_x(ids_panels{1}{1}))/(1i*omega)/dx^4;
clear Jout1_x
Jout1_y = ifftn(fJ_y .* fN_ch_all{2,2,1}+fJ_x .* conj(fN_ch_all{1,2,1})+fJ_z .* fN_ch_all{2,3,1});% yx,yy,yz
Jout1_y = Jout1_y(1:L,1:M+1,1:N);
if (num_diel > 0)
    JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = Jout1_y(ids_panels{1}{2});
else
    JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2)) = (Jout1_y(ids_panels{1}{2}))/(1i*omega)/dx^4;
end
clear Jout1_y
Jout1_z = ifftn(fJ_z .* fN_ch_all{3,3,1}+fJ_x .* conj(fN_ch_all{1,3,1})+fJ_y .* conj(fN_ch_all{2,3,1}));% zx,zy,zz
Jout1_z = Jout1_z(1:L,1:M,1:N+1);
if (num_diel > 0)
    JOut_full(N_c+inds_glob(5,1):N_c+inds_glob(5,2)) = Jout1_z(ids_panels{1}{3});
else
    JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = (Jout1_z(ids_panels{1}{3}))/(1i*omega)/dx^4;
end
clear Jout1_z

if (num_diel > 0)
    % dielectric
    Jout1_x = ifftn(fJ_x .* fN_ch_all{1,1,2} + fJ_y .* fN_ch_all{1,2,2} + fJ_z .* fN_ch_all{1,3,2});% xx,xy,xz
    Jout1_x = Jout1_x(1:L+1,1:M,1:N);
    JOut_full(N_c+inds_glob(2,1)+locs_diel_pos{1}-1) = Jout1_x(ids_panel_diel{1});
    JOut_full(N_c+inds_glob(2,1)+locs_diel_pos{4}-1) = -Jout1_x(ids_panel_diel{4});
    clear Jout1_x
    JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2))=coe_xx+JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2));
    Jout1_y = ifftn(fJ_y .* fN_ch_all{2,2,2}+fJ_x .* fN_ch_all{2,1,2}+fJ_z .* fN_ch_all{2,3,2});%yx,yy,yz
    Jout1_y = Jout1_y(1:L,1:M+1,1:N);
    JOut_full(N_c+inds_glob(4,1)+locs_diel_pos{2}-1) = Jout1_y(ids_panel_diel{2});
    JOut_full(N_c+inds_glob(4,1)+locs_diel_pos{5}-1) = -Jout1_y(ids_panel_diel{5});
    clear Jout1_y
    JOut_full(N_c+inds_glob(4,1):N_c+inds_glob(4,2))=coe_yy+JOut_full(N_c+inds_glob(4,1):N_c+inds_glob(4,2));
    Jout1_z = ifftn(fJ_z .* fN_ch_all{3,3,2}+fJ_x .* fN_ch_all{3,1,2}+fJ_y .* fN_ch_all{3,2,2});   % zx,zy,zz
    Jout1_z = Jout1_z(1:L,1:M,1:N+1);
    JOut_full(N_c+inds_glob(6,1)+locs_diel_pos{3}-1) = Jout1_z(ids_panel_diel{3});
    JOut_full(N_c+inds_glob(6,1)+locs_diel_pos{6}-1) = -Jout1_z(ids_panel_diel{6});
    clear Jout1_z
    JOut_full(N_c+inds_glob(6,1):N_c+inds_glob(6,2))=coe_zz+JOut_full(N_c+inds_glob(6,1):N_c+inds_glob(6,2));
    
end

% Jout1 = ifftn(Jout1);
% JOut(:,:,:,1) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,1) - Jout1(1:L,1:M,1:N);
% Jout2 = ifftn(Jout2);
% JOut(:,:,:,2) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,2) - Jout2(1:L,1:M,1:N);
% Jout3 = ifftn(Jout3);
% JOut(:,:,:,3) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,3) - Jout3(1:L,1:M,1:N);
% Jout4 = ifftn(Jout4);
% JOut(:,:,:,4) = (dx/6/(dx^2)) .* oneoversigma_e .* JIn(:,:,:,4) - Jout4(1:L,1:M,1:N);
% Jout5 = ifftn(Jout5);
% JOut(:,:,:,5) = (dx/2/(dx^2)) .* oneoversigma_e .* JIn(:,:,:,5) - Jout5(1:L,1:M,1:N);


% -------------------------------------------------------------------------
% Return local coordinates related to material positions
% -------------------------------------------------------------------------

if (fl_gpu == 1)
    % get from GPU
    JOut = gather(JOut(idx));
    % clear gpu data
    clear JIn; clear Jout1; clear Jout2; clear Jout3; clear fJ;
else
    JOut = JOut(idx);
end

%JOut = JOut(idx);

JOut_full(1:N_c) = JOut;

if(fl_profile == 1); disp(['Time for matvect - fft part::: ',num2str(toc)]); end;

% ---------------------------------------------------------------------
% Adding contributions due to nodal incidence matrix
% ---------------------------------------------------------------------
tic
% Perform multiplications without assigning to dum_block

JOut_full(1:N_c+N_in) = JOut_full(1:N_c+N_in) - (A_ind'*JIn0(N_c+N_in+1:N_c+2*N_in)) ;

JOut_full(N_c+N_in+1:N_c+2*N_in) = A_ind*JIn0(1:N_c+N_in);

if(fl_profile == 1); disp(['Time for matvect - A_ind matrices part::: ',num2str(toc)]); end

% ---------------------------------------------------------------------
% Sparse preconditioner [E F; G H]
% ---------------------------------------------------------------------
Ae=A_ind(:,1:N_c);
Aq=A_ind(:,N_c+1:N_c+N_in);
tic
if (fl_sparseprecon == 1)
    [JOut_full]=lse_sparse_precon_multiply(JOut_full,Ae,Aq,nodeid_4_grnd,nodeid_4_injectcurr,blk_pre,ids_pre,ord_pre);
end
if(fl_profile == 1); disp(['Time for matvect - sparse preconditioner part::: ',num2str(toc)]); end
% JOut_full(N_c+inds_glob(1,1):N_c+inds_glob(1,2)) = JOut_full(N_c+inds_glob(1,1):N_c+inds_glob(1,2))/(1i*omega)/dx^4;
% JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2)) = JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2))/(1i*omega)/dx^4;
% JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2))/(1i*omega)/dx^4;
% JOut_full_in = JOut_full;
% % block E contribution
% JOut_full(1:N_c) = A_inv2*JOut_full_in(1:N_c)+A_inv2*(-A_ind')*...
%     QQ * (UU \ (LL \ (PP * (RR \ (A_ind*A_inv2*JOut_full_in(1:N_c))))));
%
% % block F contribution
% JOut_full(1:N_c) = JOut_full(1:N_c)...
%     + A_inv2 * (A_ind') * QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(N_c+1:N_c+N_in))))));
%
%
% % block G contribution
% JOut_full(N_c+1:N_c+N_in) = ...
%     -QQ * (UU \ (LL \ (PP * (RR \ (A_ind*A_inv2*JOut_full_in(1:N_c))))));
% % block H contribution
% JOut_full(N_c+1:N_c+N_in) = JOut_full(N_c+1:N_c+N_in)...
%     +QQ * (UU \ (LL \ (PP * (RR \ (JOut_full_in(N_c+1:N_c+N_in))))));


