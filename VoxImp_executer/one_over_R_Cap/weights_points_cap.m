function [smpl_pnts,smpl_wghts] = weights_points(N,dim,intervals)
%function [W,X,Y,Z,Xp,Yp,Zp] = weights_points(N,dim,intervals)
%clc;close all;clear all;
%N=3; dim=4;  intervals =[0 1; 0 1; 0 1; 0 1];
%N=3; dim=2;  intervals =[-1 1;-1 1 ;];
%%
% dim: dimension of integration
% N:   number of points for one dimension
% returns the N^dim points according to Gauss-Legendre quadrature (X,Y,Z,X',Y',Z')
% and the associated weights W


%% Clenshaw - Curtis
% % weights
% w = zeros(N,1);
% w(1) = 2;
% w(3:2:N) = 2./(1-(2:2:N-1).^2);
% q = [w;w(N-1:-1:2)];
% w1 = real(ifft(q));
% w1d = [w1(1); 2*w1(2:N-1); w1(N)];
% 
% % points
% t = pi/(N-1)*(0:N-1)';
% x1d = cos(t);


%% Gauss-Legendre
[ w1d , x1d ] = gauss_1d_cap(N);
% points
x1d = transpose(x1d);
% weights
w1d = transpose(w1d);

int_lens=zeros(dim,1);
int_cens=zeros(dim,1);
for kk=1:dim
   int_lens(kk) = (intervals(kk,2)-intervals(kk,1))*0.5;
   int_cens(kk) = (intervals(kk,2)+intervals(kk,1))*0.5;
end


%% expansion up to 6 dimensions
if dim == 1
    
    W = w1d;
    X = x1d;
    
    smpl_pnts=[X*int_lens(1)+int_cens(1)];
    smpl_wghts=W*int_lens(1);
    
elseif dim == 2
    
    w2d = w1d*transpose(w1d);
    W = w2d(:);
    [X,Y] = meshgrid(x1d);
    X = X(:);
    Y = Y(:);
    
    smpl_pnts=[X*int_lens(1)+int_cens(1) Y*int_lens(2)+int_cens(2)];
    smpl_wghts=W*int_lens(1)*int_lens(2);
    
elseif dim == 3
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    W = w3d(:);
    [X,Y,Z] = meshgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    smpl_pnts=[X*int_lens(1)+int_cens(1) Y*int_lens(2)+int_cens(2) Z*int_lens(3)+int_cens(3)];
    smpl_wghts=W*int_lens(1)*int_lens(2)*int_lens(3);
    
elseif dim == 4
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    w4d = w3d(:)*transpose(w1d);
    W = w4d(:);
    [X,Y,Z,Xp] = ndgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Xp = Xp(:);
    
    smpl_pnts=[X*int_lens(1)+int_cens(1) Y*int_lens(2)+int_cens(2) Z*int_lens(3)+int_cens(3) Xp*int_lens(4)+int_cens(4)];
    smpl_wghts=W*int_lens(1)*int_lens(2)*int_lens(3)*int_lens(4);
    
    
elseif dim == 5
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    w4d = w3d(:)*transpose(w1d);
    w5d = w4d(:)*transpose(w1d);
    W = w5d(:);
    [X,Y,Z,Xp,Yp] = ndgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Xp = Xp(:);
    Yp = Yp(:);
    
    smpl_pnts=[X*int_lens(1)+int_cens(1) Y*int_lens(2)+int_cens(2) Z*int_lens(3)+int_cens(3) Xp*int_lens(4)+int_cens(4) Yp*int_lens(5)+int_cens(5)];
    smpl_wghts=W*int_lens(1)*int_lens(2)*int_lens(3)*int_lens(4)*int_lens(5);    
    
elseif dim == 6
    
    w2d = w1d*transpose(w1d);
    w3d = w2d(:)*transpose(w1d);
    w4d = w3d(:)*transpose(w1d);
    w5d = w4d(:)*transpose(w1d);
    w6d = w5d(:)*transpose(w1d);
    W = w6d(:);
    [X,Y,Z,Xp,Yp,Zp] = ndgrid(x1d);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Xp = Xp(:);
    Yp = Yp(:);
    Zp = Zp(:);
    
    smpl_pnts=[X*int_lens(1)+int_cens(1) Y*int_lens(2)+int_cens(2) Z*int_lens(3)+int_cens(3) Xp*int_lens(4)+int_cens(4) Yp*int_lens(5)+int_cens(5) Zp*int_lens(6)+int_cens(6)];
    smpl_wghts=W*int_lens(1)*int_lens(2)*int_lens(3)*int_lens(4)*int_lens(5)*int_lens(6);    
    
end

 
    
%end