function [JOut_full]=lse_matvect_mult1_dia(JIn0, fN_all,fN_ch_all, A_ind,dim_ind, OneoverMc, dx, freq, ko,idx,fl_sparseprecon,inds_glob,locs_diel_pos,ids_panel_diel,ids_panels,num_diel,ss,fl_Tuck)


% -------------------------------------------------------------------------
% Prepare data
% -------------------------------------------------------------------------
fJVIE=1;
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

N_c=dim_ind(1); N_n=dim_ind(2);% N_out=dim_ind(3);
% % JIn0(N_c+N_out+1:N_c+2*N_out)=JIn0(N_c+N_out+1:N_c+2*N_out)*dx^2;
% % fft dimensions
% if fl_Tuck==0
%     [LfN, MfN, NfN] = size(fN_all{1});
%     for ii=1:7
%         fN_all{ii}=fN_all{ii}*ko^2;
%     end
% else
%     LfN = fN_all{1}{5}(1); MfN = fN_all{1}{5}(2);  NfN = fN_all{1}{5}(3);
%     for ii=1:7
%         fN_all{ii}{4}=fN_all{ii}{4}*ko^2;
%     end
% end
if fl_Tuck==0
    [LfN, MfN, NfN] = size(fN_all{1});
    for ii=1:7
        if fJVIE==1
            fN_all{ii}=fN_all{ii}*ko^2;
        else
            fN_all{ii}=fN_all{ii}*jomegamu;
        end
    end
else
    LfN = fN_all{1}{5}(1); MfN = fN_all{1}{5}(2);  NfN = fN_all{1}{5}(3);
    for ii=1:7
        if fJVIE==1
            fN_all{ii}{4}=fN_all{ii}{4}*ko^2;
        else
            fN_all{ii}{4}=fN_all{ii}{4}*jomegamu;
        end
    end
end



% [LfN, MfN, NfN, ~] = size(fN_all);
if fl_Tuck==1
    LfN_x = fN_ch_all{1}(1,1); MfN_x = fN_ch_all{1}(1,2); NfN_x = fN_ch_all{1}(1,3);  %xx
    LfN_xy = fN_ch_all{1}(4,1); MfN_xy = fN_ch_all{1}(4,2); NfN_xy = fN_ch_all{1}(4,3); %xy
    LfN_xz = fN_ch_all{1}(5,1); MfN_xz = fN_ch_all{1}(5,2); NfN_xz = fN_ch_all{1}(5,3); %xz
else
    [LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
    [LfN_xy, MfN_xy, NfN_xy, ~] = size(fN_ch_all{1,2,1}); %xy
    [LfN_xz, MfN_xz, NfN_xz, ~] = size(fN_ch_all{1,3,1}); %xz
end
% [LfN_x, MfN_x, NfN_x, ~] = size(fN_ch_all{1,1,1}); %xx
% [LfN_xy, MfN_xy, NfN_xy, ~] = size(fN_ch_all{1,2,1}); %xy
% [LfN_xz, MfN_xz, NfN_xz, ~] = size(fN_ch_all{1,3,1}); %xz
% domain dimensions
[L, M, N] = size(OneoverMc);

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


JInN = zeros(L, M, N, 5);
JIn=cell(5,1);
JIn{1}=zeros(L, M, N);
JIn{2}=zeros(L, M, N);
JIn{3}=zeros(L, M, N);
JIn{4}=zeros(L, M, N);
JIn{5}=zeros(L, M, N);
JOut = zeros(L, M, N, 5);
% translate from local (idx) to global (L,M,N) coordinates
%     JInN(idx) = JIn0(1:num_curr);
JInN(idx) = JIn0(1:N_c);
JIn{1}=JInN(:,:,:,1); JIn{2}=JInN(:,:,:,2); JIn{3}=JInN(:,:,:,3); JIn{4}=JInN(:,:,:,4); JIn{5}=JInN(:,:,:,5);
clear JInN
% JOut_full = zeros(N_c+N_in+2*N_out,1);
% JOut_full_in = zeros(N_c+N_in+2*N_out,1);
JOut_full = zeros(N_c+2*N_n,1);
JOut_full_in = zeros(N_c+2*N_n,1);

% ---------------------------------------------------------------------
% apply fft and mv-op for each of the components of JIn
% ---------------------------------------------------------------------
% JOut(:,:,:,1) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,1) + Jout1(1:L,1:M,1:N);
% JOut(:,:,:,2) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,2) + Jout2(1:L,1:M,1:N);
% JOut(:,:,:,3) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,3) + Jout3(1:L,1:M,1:N);
% JOut(:,:,:,4) = (1/(6*dx)) .* oneoversigma_e .* JIn(:,:,:,4) + Jout4(1:L,1:M,1:N);
% JOut(:,:,:,5) = (1/(2*dx)) .* oneoversigma_e .* JIn(:,:,:,5) + Jout5(1:L,1:M,1:N);
if fJVIE==1
    JOut(:,:,:,1) = (1/dx) .* OneoverMc .* JIn{1};
    JOut(:,:,:,2) = (1/dx) .* OneoverMc .* JIn{2};
    JOut(:,:,:,3) = (1/dx) .* OneoverMc .* JIn{3};
    JOut(:,:,:,4) = (1/(6*dx)) .* OneoverMc .* JIn{4};
    JOut(:,:,:,5) = (1/(2*dx)) .* OneoverMc .* JIn{5};
else
    JOut(:,:,:,1) = (1/dx) .* oneoversigma_e .* JIn{1};
    JOut(:,:,:,2) = (1/dx) .* oneoversigma_e .* JIn{2};
    JOut(:,:,:,3) = (1/dx) .* oneoversigma_e .* JIn{3};
    JOut(:,:,:,4) = (1/(6*dx)) .* oneoversigma_e .* JIn{4};
    JOut(:,:,:,5) = (1/(2*dx)) .* oneoversigma_e .* JIn{5};
end

JIn{1} = fftn(JIn{1},[LfN, MfN, NfN]);
JIn{2} = fftn(JIn{2},[LfN, MfN, NfN]);
JIn{3} = fftn(JIn{3},[LfN, MfN, NfN]);
JIn{4} = fftn(JIn{4},[LfN, MfN, NfN]);
JIn{5} = fftn(JIn{5},[LfN, MfN, NfN]);
%Jout1
if fl_Tuck==0
    temp = fN_all{1} .* JIn{1}; % Gxx*Jx
    temp = temp + fN_all{2} .* (JIn{4}+JIn{5}); % Gx2d*J2d
else
    temp = ten_mat_prod(fN_all{1}{4},{fN_all{1}{1},fN_all{1}{2},fN_all{1}{3}}) .* JIn{1}; % Gxx*Jx
    temp = temp + ten_mat_prod(fN_all{2}{4},{fN_all{2}{1},fN_all{2}{2},fN_all{2}{3}}) .* (JIn{4}+JIn{5}); % Gx2d*J2d
end
% temp = temp + fN_all{2}(:,:,:) .* JIn{5}; % Gx3d*J3d
temp = ifftn(temp);
if fJVIE==0
    JOut(:,:,:,1) = JOut(:,:,:,1) + temp(1:L,1:M,1:N);
else
    JOut(:,:,:,1) = JOut(:,:,:,1) - temp(1:L,1:M,1:N);
end
%Jout2
if fl_Tuck==0
    temp = fN_all{1} .* JIn{2}; % Gyy*Jy
    temp = temp + fN_all{3} .* (JIn{4}-JIn{5}); % Gy2d*J2d - (-(y-yc))
    % temp = temp - fN_all{3}(:,:,:) .* JIn{5}; % Gy3d*J3d - (y-yc)
else
    temp = ten_mat_prod(fN_all{1}{4},{fN_all{1}{1},fN_all{1}{2},fN_all{1}{3}}) .* JIn{2}; % Gyy*Jy
    temp = temp + ten_mat_prod(fN_all{3}{4},{fN_all{3}{1},fN_all{3}{2},fN_all{3}{3}}) .* (JIn{4}-JIn{5}); % Gy2d*J2d - (-(y-yc))
    % temp = temp - fN_all{3}(:,:,:) .* JIn{5}; % Gy3d*J3d - (y-yc)
end
temp = ifftn(temp);
if fJVIE==0
    JOut(:,:,:,2) = JOut(:,:,:,2) + temp(1:L,1:M,1:N);
else
    JOut(:,:,:,2) = JOut(:,:,:,2) - temp(1:L,1:M,1:N);
end

%Jout3
if fl_Tuck==0
    temp = fN_all{1} .* JIn{3}; % Gzz*Jz
    temp = temp + fN_all{5} .* JIn{5}; % Gz3d*J3d
else
    temp = ten_mat_prod(fN_all{1}{4},{fN_all{1}{1},fN_all{1}{2},fN_all{1}{3}}) .* JIn{3}; % Gzz*Jz
    temp = temp + ten_mat_prod(fN_all{5}{4},{fN_all{5}{1},fN_all{5}{2},fN_all{5}{3}}) .* JIn{5}; % Gz3d*J3d
end
temp = ifftn(temp);
if fJVIE==0
    JOut(:,:,:,3) = JOut(:,:,:,3) + temp(1:L,1:M,1:N);
else
    JOut(:,:,:,3) = JOut(:,:,:,3) - temp(1:L,1:M,1:N);
end

%Jout4
if fl_Tuck==0
    temp = -fN_all{2} .* JIn{1}; % G2dx*Jx
    temp = temp - fN_all{3} .* JIn{2}; % G2dy*Jy - (-(y-yc))
    temp = temp + fN_all{4} .* JIn{4}; % G2d2d*J2d
    temp = temp + fN_all{6} .* JIn{5}; % G2d3d*J3d
else
    temp = -ten_mat_prod(fN_all{2}{4},{fN_all{2}{1},fN_all{2}{2},fN_all{2}{3}}) .* JIn{1}; % G2dx*Jx
    temp = temp - ten_mat_prod(fN_all{3}{4},{fN_all{3}{1},fN_all{3}{2},fN_all{3}{3}}) .* JIn{2}; % G2dy*Jy - (-(y-yc))
    temp = temp + ten_mat_prod(fN_all{4}{4},{fN_all{4}{1},fN_all{4}{2},fN_all{4}{3}}) .* JIn{4}; % G2d2d*J2d
    temp = temp + ten_mat_prod(fN_all{6}{4},{fN_all{6}{1},fN_all{6}{2},fN_all{6}{3}}) .* JIn{5}; % G2d3d*J3d
end
temp = ifftn(temp);
if fJVIE==0
    JOut(:,:,:,4) = JOut(:,:,:,4) + temp(1:L,1:M,1:N);
else
    JOut(:,:,:,4) = JOut(:,:,:,4) - temp(1:L,1:M,1:N);
end

%Jout5
if fl_Tuck==0
    temp=-fN_all{2} .* JIn{1}; % G3dx*Jx
    temp = temp + fN_all{3} .* JIn{2}; % G3dy*Jy - (y-yc)
    temp = temp - fN_all{5} .* JIn{3}; % G3dz*Jz
    temp = temp + fN_all{6} .* JIn{4}; % G3d2d*J2d
    temp = temp + fN_all{7} .* JIn{5}; % G3d3d*J3d
else
    temp=-ten_mat_prod(fN_all{2}{4},{fN_all{2}{1},fN_all{2}{2},fN_all{2}{3}}) .* JIn{1}; % G3dx*Jx
    temp = temp + ten_mat_prod(fN_all{3}{4},{fN_all{3}{1},fN_all{3}{2},fN_all{3}{3}}) .* JIn{2}; % G3dy*Jy - (y-yc)
    temp = temp - ten_mat_prod(fN_all{5}{4},{fN_all{5}{1},fN_all{5}{2},fN_all{5}{3}}) .* JIn{3}; % G3dz*Jz
    temp = temp + ten_mat_prod(fN_all{6}{4},{fN_all{6}{1},fN_all{6}{2},fN_all{6}{3}}) .* JIn{4}; % G3d2d*J2d
    temp = temp + ten_mat_prod(fN_all{7}{4},{fN_all{7}{1},fN_all{7}{2},fN_all{7}{3}}) .* JIn{5}; % G3d3d*J3d
end
temp = ifftn(temp);

if fJVIE==0
    JOut(:,:,:,5) = JOut(:,:,:,5) + temp(1:L,1:M,1:N);
else
    JOut(:,:,:,5) = JOut(:,:,:,5) - temp(1:L,1:M,1:N);
end

if fJVIE==1
    JOut(:,:,:,:) = JOut(:,:,:,:) * oneoverjomegaeo;
end
% %%%%%%%%%%%%%%%%%% Enrico code
% % x component of JIn, store contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,1),[LfN, MfN, NfN]);
% Jout1 = fN_all{1} .* fJ; % Gxx*Jx
% Jout4 = -fN_all{2} .* fJ; % G2dx*Jx
% Jout5 = Jout4; % G3dx*Jx
% 
% % y component of JIn, add contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,2),[LfN, MfN, NfN]);
% Jout2 = fN_all{1} .* fJ; % Gyy*Jy
% Jout4 = Jout4 - fN_all{3} .* fJ; % G2dy*Jy - (-(y-yc))
% Jout5 = Jout5 + fN_all{3} .* fJ; % G3dy*Jy - (y-yc)
% 
% % z component of JIn, store contribution on 2 components of Jout
% fJ = fftn(JIn(:,:,:,3),[LfN, MfN, NfN]);
% Jout3 = fN_all{1} .* fJ; % Gzz*Jz
% Jout5 = Jout5 - fN_all{5} .* fJ; % G3dz*Jz
% 
% % 2d component of JIn, add contribution on 4 components of Jout
% fJ = fftn(JIn(:,:,:,4),[LfN, MfN, NfN]);
% Jout1 = Jout1 + fN_all{2} .* fJ; % Gx2d*J2d
% Jout2 = Jout2 + fN_all{3} .* fJ; % Gy2d*J2d - (-(y-yc))
% Jout4 = Jout4 + fN_all{4} .* fJ; % G2d2d*J2d
% Jout5 = Jout5 + fN_all{6} .* fJ; % G3d2d*J2d
% 
% % 3d component of JIn, add contribution on 3 components of Jout
% fJ = fftn(JIn(:,:,:,5),[LfN, MfN, NfN]);
% Jout1 = Jout1 + fN_all{2} .* fJ; % Gx3d*J3d
% Jout2 = Jout2 - fN_all{3} .* fJ; % Gy3d*J3d - (y-yc)
% Jout3 = Jout3 + fN_all{5} .* fJ; % Gz3d*J3d
% Jout4 = Jout4 + fN_all{6} .* fJ; % G2d3d*J3d
% Jout5 = Jout5 + fN_all{7} .* fJ; % G3d3d*J3d
% 
% % apply ifft 
% Jout1 = ifftn(Jout1);
% Jout2 = ifftn(Jout2);
% Jout3 = ifftn(Jout3);
% Jout4 = ifftn(Jout4);
% Jout5 = ifftn(Jout5);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOut(:,:,:,1) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,1) + Jout1(1:L,1:M,1:N);
% JOut(:,:,:,2) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,2) + Jout2(1:L,1:M,1:N);
% JOut(:,:,:,3) = (1/dx) .* oneoversigma_e .* JIn(:,:,:,3) + Jout3(1:L,1:M,1:N);
% JOut(:,:,:,4) = (1/(6*dx)) .* oneoversigma_e .* JIn(:,:,:,4) + Jout4(1:L,1:M,1:N);
% JOut(:,:,:,5) = (1/(2*dx)) .* oneoversigma_e .* JIn(:,:,:,5) + Jout5(1:L,1:M,1:N);
% 
% % % apply ifft and add identity term
% % Jout1 = ifftn(Jout1);
% % JOut(:,:,:,1) = (1/dx) .* OneoverMc .* JIn(:,:,:,1) - Jout1(1:L,1:M,1:N);
% % Jout2 = ifftn(Jout2);
% % JOut(:,:,:,2) = (1/dx) .* OneoverMc .* JIn(:,:,:,2) - Jout2(1:L,1:M,1:N);
% % Jout3 = ifftn(Jout3);
% % JOut(:,:,:,3) = (1/dx) .* OneoverMc .* JIn(:,:,:,3) - Jout3(1:L,1:M,1:N);
% % Jout4 = ifftn(Jout4);
% % JOut(:,:,:,4) = (dx/6/(dx^2)) .* OneoverMc .* JIn(:,:,:,4) - Jout4(1:L,1:M,1:N);
% % Jout5 = ifftn(Jout5);
% % JOut(:,:,:,5) = (dx/2/(dx^2)) .* OneoverMc .* JIn(:,:,:,5) - Jout5(1:L,1:M,1:N);
% % % multiply by 1/(jweps0)
% % JOut(:,:,:,:) = JOut(:,:,:,:) * oneoverjomegaeo;


% conductor
if fl_Tuck==1
    Jout1_x = ifftn(fJ_x .* ten_mat_prod(fN_ch_all{2},{fN_ch_all{3},fN_ch_all{4},fN_ch_all{5}}) + ...
            fJ_y .* ten_mat_prod(fN_ch_all{14},{fN_ch_all{15},fN_ch_all{16},fN_ch_all{17}}) + ...
            fJ_z .* ten_mat_prod(fN_ch_all{18},{fN_ch_all{19},fN_ch_all{20},fN_ch_all{21}})); % xx,xy,xz
    Jout1_x = Jout1_x(1:L+1,1:M,1:N);
    JOut_full(N_c+inds_glob(1,1):N_c+inds_glob(1,2)) = Jout1_x(ids_panels{1}{1})/(1i*omega)/dx^4;
    clear Jout1_x
    
    Jout1_y = ifftn(fJ_x .* ten_mat_prod(fN_ch_all{26},{fN_ch_all{27},fN_ch_all{28},fN_ch_all{29}}) + ...
            fJ_y .* ten_mat_prod(fN_ch_all{6},{fN_ch_all{7},fN_ch_all{8},fN_ch_all{9}}) + ...
            fJ_z .* ten_mat_prod(fN_ch_all{22},{fN_ch_all{23},fN_ch_all{24},fN_ch_all{25}})); % yx,yy,yz
    Jout1_y = Jout1_y(1:L,1:M+1,1:N);
    if (num_diel > 0)
        JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = Jout1_y(ids_panels{1}{2});
    else
        JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2)) = Jout1_y(ids_panels{1}{2})/(1i*omega)/dx^4;
    end
    clear Jout1_y
    
    Jout1_z = ifftn(fJ_x .* ten_mat_prod(fN_ch_all{30},{fN_ch_all{31},fN_ch_all{32},fN_ch_all{33}}) + ...
            fJ_y .* ten_mat_prod(fN_ch_all{34},{fN_ch_all{35},fN_ch_all{36},fN_ch_all{37}}) + ...
            fJ_z .* ten_mat_prod(fN_ch_all{10},{fN_ch_all{11},fN_ch_all{12},fN_ch_all{13}})); % zx,zy,zz
    Jout1_z = Jout1_z(1:L,1:M,1:N+1);
    if (num_diel > 0)
        JOut_full(N_c+inds_glob(5,1):N_c+inds_glob(5,2)) = Jout1_z(ids_panels{1}{3});
    else
        JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = Jout1_z(ids_panels{1}{3})/(1i*omega)/dx^4;
    end
    
    clear Jout1_z 
       
else
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
end



% Jout1_x = ifftn(fJ_x .* fN_ch_all{1,1,1}+fJ_y .* fN_ch_all{1,2,1}+fJ_z .* fN_ch_all{1,3,1});% xx,xy,xz
% Jout1_x = Jout1_x(1:L+1,1:M,1:N);
% JOut_full(N_c+inds_glob(1,1):N_c+inds_glob(1,2)) = (Jout1_x(ids_panels{1}{1}))/(1i*omega)/dx^4;
% clear Jout1_x
% Jout1_y = ifftn(fJ_y .* fN_ch_all{2,2,1}+fJ_x .* conj(fN_ch_all{1,2,1})+fJ_z .* fN_ch_all{2,3,1});% yx,yy,yz
% Jout1_y = Jout1_y(1:L,1:M+1,1:N);
% if (num_diel > 0)
%     JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = Jout1_y(ids_panels{1}{2});
% else
%     JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2)) = (Jout1_y(ids_panels{1}{2}))/(1i*omega)/dx^4;
% end
% clear Jout1_y
% Jout1_z = ifftn(fJ_z .* fN_ch_all{3,3,1}+fJ_x .* conj(fN_ch_all{1,3,1})+fJ_y .* conj(fN_ch_all{2,3,1}));% zx,zy,zz
% Jout1_z = Jout1_z(1:L,1:M,1:N+1);
% if (num_diel > 0)
%     JOut_full(N_c+inds_glob(5,1):N_c+inds_glob(5,2)) = Jout1_z(ids_panels{1}{3});
% else
%     JOut_full(N_c+inds_glob(3,1):N_c+inds_glob(3,2)) = (Jout1_z(ids_panels{1}{3}))/(1i*omega)/dx^4;
% end
% clear Jout1_z

if (num_diel > 0)
    % dielectric
    if fl_Tuck==1
        Jout1_x = ifftn(fJ_x .* ten_mat_prod(fN_ch_all{38},{fN_ch_all{39},fN_ch_all{40},fN_ch_all{41}}) + ...
            fJ_y .* ten_mat_prod(fN_ch_all{50},{fN_ch_all{51},fN_ch_all{52},fN_ch_all{53}}) + ...
            fJ_z .* ten_mat_prod(fN_ch_all{54},{fN_ch_all{55},fN_ch_all{56},fN_ch_all{57}})); % xx,xy,xz
        Jout1_x = Jout1_x(1:L+1,1:M,1:N);
        JOut_full(N_c+inds_glob(2,1)+locs_diel_pos{1}-1) = Jout1_x(ids_panel_diel{1});
        JOut_full(N_c+inds_glob(2,1)+locs_diel_pos{4}-1) = -Jout1_x(ids_panel_diel{4});
        clear Jout1_x
        JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2))=coe_xx+JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2));
        Jout1_y = ifftn(fJ_y .* ten_mat_prod(fN_ch_all{42},{fN_ch_all{43},fN_ch_all{44},fN_ch_all{45}})+ ...
            fJ_x .* ten_mat_prod(fN_ch_all{62},{fN_ch_all{63},fN_ch_all{64},fN_ch_all{65}})+ ...
            fJ_z .* ten_mat_prod(fN_ch_all{58},{fN_ch_all{59},fN_ch_all{60},fN_ch_all{61}})); %yx,yy,yz
        Jout1_y = Jout1_y(1:L,1:M+1,1:N);
        JOut_full(N_c+inds_glob(4,1)+locs_diel_pos{2}-1) = Jout1_y(ids_panel_diel{2});
        JOut_full(N_c+inds_glob(4,1)+locs_diel_pos{5}-1) = -Jout1_y(ids_panel_diel{5});
        clear Jout1_y
        JOut_full(N_c+inds_glob(4,1):N_c+inds_glob(4,2))=coe_yy+JOut_full(N_c+inds_glob(4,1):N_c+inds_glob(4,2));
        Jout1_z = ifftn(fJ_z .* ten_mat_prod(fN_ch_all{46},{fN_ch_all{47},fN_ch_all{48},fN_ch_all{49}})+ ...
            fJ_x .* ten_mat_prod(fN_ch_all{66},{fN_ch_all{67},fN_ch_all{68},fN_ch_all{69}})+ ...
            fJ_y .* ten_mat_prod(fN_ch_all{70},{fN_ch_all{71},fN_ch_all{72},fN_ch_all{73}}));  % zx,zy,zz
        Jout1_z = Jout1_z(1:L,1:M,1:N+1);
        JOut_full(N_c+inds_glob(6,1)+locs_diel_pos{3}-1) = Jout1_z(ids_panel_diel{3});
        JOut_full(N_c+inds_glob(6,1)+locs_diel_pos{6}-1) = -Jout1_z(ids_panel_diel{6});
        clear Jout1_z
        JOut_full(N_c+inds_glob(6,1):N_c+inds_glob(6,2))=coe_zz+JOut_full(N_c+inds_glob(6,1):N_c+inds_glob(6,2));
    else
        
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
%     Jout1_x = ifftn(fJ_x .* fN_ch_all{1,1,2} + fJ_y .* fN_ch_all{1,2,2} + fJ_z .* fN_ch_all{1,3,2});% xx,xy,xz
%     Jout1_x = Jout1_x(1:L+1,1:M,1:N);
%     JOut_full(N_c+inds_glob(2,1)+locs_diel_pos{1}-1) = Jout1_x(ids_panel_diel{1});
%     JOut_full(N_c+inds_glob(2,1)+locs_diel_pos{4}-1) = -Jout1_x(ids_panel_diel{4});
%     clear Jout1_x
%     JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2))=coe_xx+JOut_full(N_c+inds_glob(2,1):N_c+inds_glob(2,2));
%     Jout1_y = ifftn(fJ_y .* fN_ch_all{2,2,2}+fJ_x .* fN_ch_all{2,1,2}+fJ_z .* fN_ch_all{2,3,2});%yx,yy,yz
%     Jout1_y = Jout1_y(1:L,1:M+1,1:N);
%     JOut_full(N_c+inds_glob(4,1)+locs_diel_pos{2}-1) = Jout1_y(ids_panel_diel{2});
%     JOut_full(N_c+inds_glob(4,1)+locs_diel_pos{5}-1) = -Jout1_y(ids_panel_diel{5});
%     clear Jout1_y
%     JOut_full(N_c+inds_glob(4,1):N_c+inds_glob(4,2))=coe_yy+JOut_full(N_c+inds_glob(4,1):N_c+inds_glob(4,2));
%     Jout1_z = ifftn(fJ_z .* fN_ch_all{3,3,2}+fJ_x .* fN_ch_all{3,1,2}+fJ_y .* fN_ch_all{3,2,2});   % zx,zy,zz
%     Jout1_z = Jout1_z(1:L,1:M,1:N+1);
%     JOut_full(N_c+inds_glob(6,1)+locs_diel_pos{3}-1) = Jout1_z(ids_panel_diel{3});
%     JOut_full(N_c+inds_glob(6,1)+locs_diel_pos{6}-1) = -Jout1_z(ids_panel_diel{6});
%     clear Jout1_z
%     JOut_full(N_c+inds_glob(6,1):N_c+inds_glob(6,2))=coe_zz+JOut_full(N_c+inds_glob(6,1):N_c+inds_glob(6,2));
    
end


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

JOut_full(1:N_c+N_n) = JOut_full(1:N_c+N_n) - (A_ind'*JIn0(N_c+N_n+1:N_c+2*N_n)) ;

JOut_full(N_c+N_n+1:N_c+2*N_n) = A_ind*JIn0(1:N_c+N_n);

if(fl_profile == 1); disp(['Time for matvect - A_ind matrices part::: ',num2str(toc)]); end

% ---------------------------------------------------------------------
% Sparse preconditioner [E F; G H]
% ---------------------------------------------------------------------
Ae=A_ind(:,1:N_c);
Aq=A_ind(:,N_c+1:N_c+N_n);
tic
if (fl_sparseprecon == 1)
    [JOut_full]=lse_sparse_precon_multiply_dia(JOut_full,Ae,Aq);
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


