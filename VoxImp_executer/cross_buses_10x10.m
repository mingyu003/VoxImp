
Res = 0.25e-6; % voxel size (deltax)
% -------------------------------------------------------------------------
%                  Inputs for the Structure
% -------------------------------------------------------------------------
% We only need centers (Cnt), dimensions (Dims), and orientations (Orients)
% of the conductors at the end of this part.

% inputs for generating conductors with specified lengths and widths of arms
num_conds = 100; % number of conductors
num_diels = 0;
num_ports = 100; % number of ports

len=25e-6; wid=1e-6; hei=1e-6; dist=2e-6;
Cnt=zeros(100,5);
for ii=1:10
    Cnt(ii,:)=[3e-6+wid/2+dist*(ii-1) len/2 hei/2 ii 0];
    Cnt(ii+10,:)=[len/2 3e-6+wid/2+dist*(ii-1) hei/2+dist ii+10 0];
    Cnt(ii+2*10,:)=[3e-6+wid/2+dist*(ii-1) len/2 hei/2+2*dist ii+20 0];
    Cnt(ii+3*10,:)=[len/2 3e-6+wid/2+dist*(ii-1) hei/2+3*dist ii+30 0];
    Cnt(ii+4*10,:)=[3e-6+wid/2+dist*(ii-1) len/2 hei/2+4*dist ii+40 0];
    Cnt(ii+5*10,:)=[len/2 3e-6+wid/2+dist*(ii-1) hei/2+5*dist ii+50 0];
    Cnt(ii+6*10,:)=[3e-6+wid/2+dist*(ii-1) len/2 hei/2+6*dist ii+60 0];
    Cnt(ii+7*10,:)=[len/2 3e-6+wid/2+dist*(ii-1) hei/2+7*dist ii+70 0];
    Cnt(ii+8*10,:)=[3e-6+wid/2+dist*(ii-1) len/2 hei/2+8*dist ii+80 0];
    Cnt(ii+9*10,:)=[len/2 3e-6+wid/2+dist*(ii-1) hei/2+9*dist ii+90 0];
end

Dims=zeros(100,3);
for ii=1:100
    Dims(ii,:)=[len,wid,hei];
end
Orients=[];
for ii=1:10
    Orients=[Orients;'y'];
end
for ii=11:20
    Orients=[Orients;'x'];
end
Orients=[Orients;Orients;Orients;Orients;Orients];

Eps_inout=zeros(num_conds+num_diels,2);% inner and outer eps
for ii=1:num_conds
    Eps_inout(ii,:)=[0 1];
end
% -------------------------------------------------------------------------
%                  Input for Computational Domain
% -------------------------------------------------------------------------
% At the end of this part, we only need bbox_min(3) and bbox_max(3) vectors
% define computational domain or bounding box enclosing the structure
bbox_min=[0 0 0]; % minimum coordinates of bounding box (bbox) - set to positive reals if possible
bbox_max=[len len 19*hei]+1e-12; % max coordinates of bbox

% -------------------------------------------------------------------------
%                  Input for Ports
% -------------------------------------------------------------------------
% % At the end of this part, we need structures pnt_lft{xx} and pnt_rght{xx} which
% % contains the coordinates of nodes on both sides of xxth port
% port_direction=zeros(30,2);%normal direction of left and right port. 1-x 2-y 3-z
% for ii=1:30
%     port_direction(ii,:)=[1 1];
% end
% defining the nodes in first port
pnt_lft=cell(num_ports,1); %exc
pnt_rght=cell(num_ports,1); %gnd
pnt_lft{1}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{1}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{1}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 0 (2*ll-1)*(0.5*Res) 2]; %1-x 2-y 3-z
        pnt_lft{1}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 len (2*ll-1)*(0.5*Res) 2]; %1-x 2-y 3-z
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk}(:,1)=pnt_lft{1}(:,1) + dist*(kk-1);
    pnt_lft{kk}(:,[2 3 4])=pnt_lft{1}(:,[2 3 4]);
    pnt_rght{kk}(:,1)=pnt_rght{1}(:,1) + dist*(kk-1);
    pnt_rght{kk}(:,[2 3 4])=pnt_rght{1}(:,[2 3 4]);
end

pnt_lft{11}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{11}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{11}(dum,1:4)=[0 (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+dist 1]; 
        pnt_lft{11}(dum,1:4)=[len (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+dist 1]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+10}(:,2)=pnt_lft{11}(:,2) + dist*(kk-1);
    pnt_lft{kk+10}(:,[1 3 4])=pnt_lft{11}(:,[1 3 4]);
    pnt_rght{kk+10}(:,2)=pnt_rght{11}(:,2) + dist*(kk-1);
    pnt_rght{kk+10}(:,[1 3 4])=pnt_rght{11}(:,[1 3 4]);
end

pnt_lft{21}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{21}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{21}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 0 (2*ll-1)*(0.5*Res)+2*dist 2]; 
        pnt_lft{21}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 len (2*ll-1)*(0.5*Res)+2*dist 2]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+20}(:,1)=pnt_lft{21}(:,1) + dist*(kk-1);
    pnt_lft{kk+20}(:,[2 3 4])=pnt_lft{21}(:,[2 3 4]);
    pnt_rght{kk+20}(:,1)=pnt_rght{21}(:,1) + dist*(kk-1);
    pnt_rght{kk+20}(:,[2 3 4])=pnt_rght{21}(:,[2 3 4]);
end

pnt_lft{31}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{31}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{31}(dum,1:4)=[0 (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+3*dist 1]; 
        pnt_lft{31}(dum,1:4)=[len (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+3*dist 1]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+30}(:,2)=pnt_lft{31}(:,2) + dist*(kk-1);
    pnt_lft{kk+30}(:,[1 3 4])=pnt_lft{31}(:,[1 3 4]);
    pnt_rght{kk+30}(:,2)=pnt_rght{31}(:,2) + dist*(kk-1);
    pnt_rght{kk+30}(:,[1 3 4])=pnt_rght{31}(:,[1 3 4]);
end

pnt_lft{41}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{41}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{41}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 0 (2*ll-1)*(0.5*Res)+4*dist 2]; 
        pnt_lft{41}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 len (2*ll-1)*(0.5*Res)+4*dist 2]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+40}(:,1)=pnt_lft{41}(:,1) + dist*(kk-1);
    pnt_lft{kk+40}(:,[2 3 4])=pnt_lft{41}(:,[2 3 4]);
    pnt_rght{kk+40}(:,1)=pnt_rght{41}(:,1) + dist*(kk-1);
    pnt_rght{kk+40}(:,[2 3 4])=pnt_rght{41}(:,[2 3 4]);
end

pnt_lft{51}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{51}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{51}(dum,1:4)=[0 (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+5*dist 1]; 
        pnt_lft{51}(dum,1:4)=[len (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+5*dist 1]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+50}(:,2)=pnt_lft{51}(:,2) + dist*(kk-1);
    pnt_lft{kk+50}(:,[1 3 4])=pnt_lft{51}(:,[1 3 4]);
    pnt_rght{kk+50}(:,2)=pnt_rght{51}(:,2) + dist*(kk-1);
    pnt_rght{kk+50}(:,[1 3 4])=pnt_rght{51}(:,[1 3 4]);
end

pnt_lft{61}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{61}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{61}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 0 (2*ll-1)*(0.5*Res)+6*dist 2]; 
        pnt_lft{61}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 len (2*ll-1)*(0.5*Res)+6*dist 2]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+60}(:,1)=pnt_lft{61}(:,1) + dist*(kk-1);
    pnt_lft{kk+60}(:,[2 3 4])=pnt_lft{61}(:,[2 3 4]);
    pnt_rght{kk+60}(:,1)=pnt_rght{61}(:,1) + dist*(kk-1);
    pnt_rght{kk+60}(:,[2 3 4])=pnt_rght{61}(:,[2 3 4]);
end

pnt_lft{71}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{71}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{71}(dum,1:4)=[0 (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+7*dist 1]; 
        pnt_lft{71}(dum,1:4)=[len (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+7*dist 1]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+70}(:,2)=pnt_lft{71}(:,2) + dist*(kk-1);
    pnt_lft{kk+70}(:,[1 3 4])=pnt_lft{71}(:,[1 3 4]);
    pnt_rght{kk+70}(:,2)=pnt_rght{71}(:,2) + dist*(kk-1);
    pnt_rght{kk+70}(:,[1 3 4])=pnt_rght{71}(:,[1 3 4]);
end

pnt_lft{81}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{81}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{81}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 0 (2*ll-1)*(0.5*Res)+8*dist 2]; 
        pnt_lft{81}(dum,1:4)=[(2*kk-1)*(0.5*Res)+3e-6 len (2*ll-1)*(0.5*Res)+8*dist 2]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+80}(:,1)=pnt_lft{81}(:,1) + dist*(kk-1);
    pnt_lft{kk+80}(:,[2 3 4])=pnt_lft{81}(:,[2 3 4]);
    pnt_rght{kk+80}(:,1)=pnt_rght{81}(:,1) + dist*(kk-1);
    pnt_rght{kk+80}(:,[2 3 4])=pnt_rght{81}(:,[2 3 4]);
end

pnt_lft{91}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{91}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;
for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{91}(dum,1:4)=[0 (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+9*dist 1]; 
        pnt_lft{91}(dum,1:4)=[len (2*kk-1)*(0.5*Res)+3e-6 (2*ll-1)*(0.5*Res)+9*dist 1]; 
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:10
    pnt_lft{kk+90}(:,2)=pnt_lft{91}(:,2) + dist*(kk-1);
    pnt_lft{kk+90}(:,[1 3 4])=pnt_lft{91}(:,[1 3 4]);
    pnt_rght{kk+90}(:,2)=pnt_rght{91}(:,2) + dist*(kk-1);
    pnt_rght{kk+90}(:,[1 3 4])=pnt_rght{91}(:,[1 3 4]);
end

% defining nodes connected ground if conductors without ports exist; if
% there is no, then leave as a empty array.

%pnt_well_cond=[pnt_lft{1}(1,1) pnt_lft{1}(1,2) + dist pnt_lft{1}(1,3);];
pnt_well_cond=[];

