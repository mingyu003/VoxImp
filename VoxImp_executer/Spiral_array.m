
% % -------------------------------------------------------------------------
% %                  Inputs for Simulation
% % -------------------------------------------------------------------------
Res = 0.2e-6; % voxel size (deltax)
% -------------------------------------------------------------------------
%                  Inputs for the Structure
% -------------------------------------------------------------------------
% We only need centers (Cnt), dimensions (Dims), and orientations (Orients)
% of the conductors at the end of this part.

% inputs for generating conductors with specified lengths and widths of arms
num_conds = 128; % number of conductors
num_diels = 0;
num_ports = 128; % number of ports

ilen=12e-6; % the length of initial arm
wid=1.2e-6; hei=1.2e-6; 
spa=0.6e-6;% space of two arms


len1=ilen+5*spa+5*wid; len2=len1; len3=ilen+4*spa+4*wid; len4=len3; len5=ilen+3*spa+3*wid; len6=len5;
len7=ilen+2*spa+2*wid; len8=len7; len9=ilen+spa+wid; len10=len9; len11=ilen; len12=len11;

gap=2.4e-6; % distance btw adjacent two spirals
gapx=gap+len1; gapy=gap+wid+len2;
% cond 1-16
cen1_1=[len1/2 wid/2 hei/2 1 0]; %4th col:cond numbering, 0 for GND or diel, %5th col: diel numbering, 0 for no diel
cen1_5=[spa+wid+len5/2 wid+spa+wid/2 hei/2 1 0]; cen1_9=[2*spa+2*wid+len9/2 2*wid+2*spa+wid/2 hei/2 1 0];
cen1_2=[len1-wid/2 wid+len2/2 hei/2 1 0]; cen1_6=[spa+wid+len5-wid/2 2*wid+spa+len6/2 hei/2 1 0]; cen1_10=[len5-wid/2 3*wid+2*spa+len10/2 hei/2 1 0];
cen1_3=[len1-wid-len3/2 wid+len2-wid/2 hei/2 1 0]; cen1_7=[2*spa+wid+len7/2 2*wid+spa+len6-wid/2 hei/2 1 0]; cen1_11=[3*spa+2*wid+len11/2 wid+len6-wid/2 hei/2 1 0];
cen1_4=[spa+wid/2 wid+spa+len4/2 hei/2 1 0]; cen1_8=[2*spa+wid+wid/2 2*wid+2*spa+len8/2 hei/2 1 0]; cen1_12=[3*spa+2*wid+wid/2 3*wid+3*spa+len12/2 hei/2 1 0];
cnt_cond1=[cen1_1; cen1_2; cen1_3; cen1_4; cen1_5; cen1_6; cen1_7; cen1_8; cen1_9; cen1_10; cen1_11; cen1_12];
cnt_cond2=cnt_cond1; cnt_cond2(:,1)=cnt_cond2(:,1)+gapx; cnt_cond2(:,4)=2;
cnt_cond3=cnt_cond1; cnt_cond3(:,1)=cnt_cond3(:,1)+2*gapx; cnt_cond3(:,4)=3;
cnt_cond4=cnt_cond1; cnt_cond4(:,1)=cnt_cond4(:,1)+3*gapx; cnt_cond4(:,4)=4;
cnt_cond5=cnt_cond1; cnt_cond5(:,1)=cnt_cond5(:,1)+4*gapx; cnt_cond5(:,4)=5;
cnt_cond6=cnt_cond1; cnt_cond6(:,1)=cnt_cond6(:,1)+5*gapx; cnt_cond6(:,4)=6;
cnt_cond7=cnt_cond1; cnt_cond7(:,1)=cnt_cond7(:,1)+6*gapx; cnt_cond7(:,4)=7;
cnt_cond8=cnt_cond1; cnt_cond8(:,1)=cnt_cond8(:,1)+7*gapx; cnt_cond8(:,4)=8;
cnt_cond9=cnt_cond1; cnt_cond9(:,1)=cnt_cond9(:,1)+8*gapx; cnt_cond9(:,4)=9;
cnt_cond10=cnt_cond1; cnt_cond10(:,1)=cnt_cond10(:,1)+9*gapx; cnt_cond10(:,4)=10;
cnt_cond11=cnt_cond1; cnt_cond11(:,1)=cnt_cond11(:,1)+10*gapx; cnt_cond11(:,4)=11;
cnt_cond12=cnt_cond1; cnt_cond12(:,1)=cnt_cond12(:,1)+11*gapx; cnt_cond12(:,4)=12;
cnt_cond13=cnt_cond1; cnt_cond13(:,1)=cnt_cond13(:,1)+12*gapx; cnt_cond13(:,4)=13;
cnt_cond14=cnt_cond1; cnt_cond14(:,1)=cnt_cond14(:,1)+13*gapx; cnt_cond14(:,4)=14;
cnt_cond15=cnt_cond1; cnt_cond15(:,1)=cnt_cond15(:,1)+14*gapx; cnt_cond15(:,4)=15;
cnt_cond16=cnt_cond1; cnt_cond16(:,1)=cnt_cond16(:,1)+15*gapx; cnt_cond16(:,4)=16;
% cond 17-32
cnt_cond17=cnt_cond1; cnt_cond17(:,2)=cnt_cond17(:,2)+gapy; cnt_cond17(:,4)=17;
cnt_cond18=cnt_cond2; cnt_cond18(:,2)=cnt_cond18(:,2)+gapy; cnt_cond18(:,4)=18;
cnt_cond19=cnt_cond3; cnt_cond19(:,2)=cnt_cond19(:,2)+gapy; cnt_cond19(:,4)=19;
cnt_cond20=cnt_cond4; cnt_cond20(:,2)=cnt_cond20(:,2)+gapy; cnt_cond20(:,4)=20;
cnt_cond21=cnt_cond5; cnt_cond21(:,2)=cnt_cond21(:,2)+gapy; cnt_cond21(:,4)=21;
cnt_cond22=cnt_cond6; cnt_cond22(:,2)=cnt_cond22(:,2)+gapy; cnt_cond22(:,4)=22;
cnt_cond23=cnt_cond7; cnt_cond23(:,2)=cnt_cond23(:,2)+gapy; cnt_cond23(:,4)=23;
cnt_cond24=cnt_cond8; cnt_cond24(:,2)=cnt_cond24(:,2)+gapy; cnt_cond24(:,4)=24;
cnt_cond25=cnt_cond9; cnt_cond25(:,2)=cnt_cond25(:,2)+gapy; cnt_cond25(:,4)=25;
cnt_cond26=cnt_cond10; cnt_cond26(:,2)=cnt_cond26(:,2)+gapy; cnt_cond26(:,4)=26;
cnt_cond27=cnt_cond11; cnt_cond27(:,2)=cnt_cond27(:,2)+gapy; cnt_cond27(:,4)=27;
cnt_cond28=cnt_cond12; cnt_cond28(:,2)=cnt_cond28(:,2)+gapy; cnt_cond28(:,4)=28;
cnt_cond29=cnt_cond13; cnt_cond29(:,2)=cnt_cond29(:,2)+gapy; cnt_cond29(:,4)=29;
cnt_cond30=cnt_cond14; cnt_cond30(:,2)=cnt_cond30(:,2)+gapy; cnt_cond30(:,4)=30;
cnt_cond31=cnt_cond15; cnt_cond31(:,2)=cnt_cond31(:,2)+gapy; cnt_cond31(:,4)=31;
cnt_cond32=cnt_cond16; cnt_cond32(:,2)=cnt_cond32(:,2)+gapy; cnt_cond32(:,4)=32;
% cond 33-48
cnt_cond33=cnt_cond1; cnt_cond33(:,2)=cnt_cond33(:,2)+2*gapy; cnt_cond33(:,4)=33;
cnt_cond34=cnt_cond2; cnt_cond34(:,2)=cnt_cond34(:,2)+2*gapy; cnt_cond34(:,4)=34;
cnt_cond35=cnt_cond3; cnt_cond35(:,2)=cnt_cond35(:,2)+2*gapy; cnt_cond35(:,4)=35;
cnt_cond36=cnt_cond4; cnt_cond36(:,2)=cnt_cond36(:,2)+2*gapy; cnt_cond36(:,4)=36;
cnt_cond37=cnt_cond5; cnt_cond37(:,2)=cnt_cond37(:,2)+2*gapy; cnt_cond37(:,4)=37;
cnt_cond38=cnt_cond6; cnt_cond38(:,2)=cnt_cond38(:,2)+2*gapy; cnt_cond38(:,4)=38;
cnt_cond39=cnt_cond7; cnt_cond39(:,2)=cnt_cond39(:,2)+2*gapy; cnt_cond39(:,4)=39;
cnt_cond40=cnt_cond8; cnt_cond40(:,2)=cnt_cond40(:,2)+2*gapy; cnt_cond40(:,4)=40;
cnt_cond41=cnt_cond9; cnt_cond41(:,2)=cnt_cond41(:,2)+2*gapy; cnt_cond41(:,4)=41;
cnt_cond42=cnt_cond10; cnt_cond42(:,2)=cnt_cond42(:,2)+2*gapy; cnt_cond42(:,4)=42;
cnt_cond43=cnt_cond11; cnt_cond43(:,2)=cnt_cond43(:,2)+2*gapy; cnt_cond43(:,4)=43;
cnt_cond44=cnt_cond12; cnt_cond44(:,2)=cnt_cond44(:,2)+2*gapy; cnt_cond44(:,4)=44;
cnt_cond45=cnt_cond13; cnt_cond45(:,2)=cnt_cond45(:,2)+2*gapy; cnt_cond45(:,4)=45;
cnt_cond46=cnt_cond14; cnt_cond46(:,2)=cnt_cond46(:,2)+2*gapy; cnt_cond46(:,4)=46;
cnt_cond47=cnt_cond15; cnt_cond47(:,2)=cnt_cond47(:,2)+2*gapy; cnt_cond47(:,4)=47;
cnt_cond48=cnt_cond16; cnt_cond48(:,2)=cnt_cond48(:,2)+2*gapy; cnt_cond48(:,4)=48;
% cond 49-64
cnt_cond49=cnt_cond1; cnt_cond49(:,2)=cnt_cond49(:,2)+3*gapy; cnt_cond49(:,4)=49;
cnt_cond50=cnt_cond2; cnt_cond50(:,2)=cnt_cond50(:,2)+3*gapy; cnt_cond50(:,4)=50;
cnt_cond51=cnt_cond3; cnt_cond51(:,2)=cnt_cond51(:,2)+3*gapy; cnt_cond51(:,4)=51;
cnt_cond52=cnt_cond4; cnt_cond52(:,2)=cnt_cond52(:,2)+3*gapy; cnt_cond52(:,4)=52;
cnt_cond53=cnt_cond5; cnt_cond53(:,2)=cnt_cond53(:,2)+3*gapy; cnt_cond53(:,4)=53;
cnt_cond54=cnt_cond6; cnt_cond54(:,2)=cnt_cond54(:,2)+3*gapy; cnt_cond54(:,4)=54;
cnt_cond55=cnt_cond7; cnt_cond55(:,2)=cnt_cond55(:,2)+3*gapy; cnt_cond55(:,4)=55;
cnt_cond56=cnt_cond8; cnt_cond56(:,2)=cnt_cond56(:,2)+3*gapy; cnt_cond56(:,4)=56;
cnt_cond57=cnt_cond9; cnt_cond57(:,2)=cnt_cond57(:,2)+3*gapy; cnt_cond57(:,4)=57;
cnt_cond58=cnt_cond10; cnt_cond58(:,2)=cnt_cond58(:,2)+3*gapy; cnt_cond58(:,4)=58;
cnt_cond59=cnt_cond11; cnt_cond59(:,2)=cnt_cond59(:,2)+3*gapy; cnt_cond59(:,4)=59;
cnt_cond60=cnt_cond12; cnt_cond60(:,2)=cnt_cond60(:,2)+3*gapy; cnt_cond60(:,4)=60;
cnt_cond61=cnt_cond13; cnt_cond61(:,2)=cnt_cond61(:,2)+3*gapy; cnt_cond61(:,4)=61;
cnt_cond62=cnt_cond14; cnt_cond62(:,2)=cnt_cond62(:,2)+3*gapy; cnt_cond62(:,4)=62;
cnt_cond63=cnt_cond15; cnt_cond63(:,2)=cnt_cond63(:,2)+3*gapy; cnt_cond63(:,4)=63;
cnt_cond64=cnt_cond16; cnt_cond64(:,2)=cnt_cond64(:,2)+3*gapy; cnt_cond64(:,4)=64;
% cond 65-80
cnt_cond65=cnt_cond1; cnt_cond65(:,2)=cnt_cond65(:,2)+4*gapy; cnt_cond65(:,4)=65;
cnt_cond66=cnt_cond2; cnt_cond66(:,2)=cnt_cond66(:,2)+4*gapy; cnt_cond66(:,4)=66;
cnt_cond67=cnt_cond3; cnt_cond67(:,2)=cnt_cond67(:,2)+4*gapy; cnt_cond67(:,4)=67;
cnt_cond68=cnt_cond4; cnt_cond68(:,2)=cnt_cond68(:,2)+4*gapy; cnt_cond68(:,4)=68;
cnt_cond69=cnt_cond5; cnt_cond69(:,2)=cnt_cond69(:,2)+4*gapy; cnt_cond69(:,4)=69;
cnt_cond70=cnt_cond6; cnt_cond70(:,2)=cnt_cond70(:,2)+4*gapy; cnt_cond70(:,4)=70;
cnt_cond71=cnt_cond7; cnt_cond71(:,2)=cnt_cond71(:,2)+4*gapy; cnt_cond71(:,4)=71;
cnt_cond72=cnt_cond8; cnt_cond72(:,2)=cnt_cond72(:,2)+4*gapy; cnt_cond72(:,4)=72;
cnt_cond73=cnt_cond9; cnt_cond73(:,2)=cnt_cond73(:,2)+4*gapy; cnt_cond73(:,4)=73;
cnt_cond74=cnt_cond10; cnt_cond74(:,2)=cnt_cond74(:,2)+4*gapy; cnt_cond74(:,4)=74;
cnt_cond75=cnt_cond11; cnt_cond75(:,2)=cnt_cond75(:,2)+4*gapy; cnt_cond75(:,4)=75;
cnt_cond76=cnt_cond12; cnt_cond76(:,2)=cnt_cond76(:,2)+4*gapy; cnt_cond76(:,4)=76;
cnt_cond77=cnt_cond13; cnt_cond77(:,2)=cnt_cond77(:,2)+4*gapy; cnt_cond77(:,4)=77;
cnt_cond78=cnt_cond14; cnt_cond78(:,2)=cnt_cond78(:,2)+4*gapy; cnt_cond78(:,4)=78;
cnt_cond79=cnt_cond15; cnt_cond79(:,2)=cnt_cond79(:,2)+4*gapy; cnt_cond79(:,4)=79;
cnt_cond80=cnt_cond16; cnt_cond80(:,2)=cnt_cond80(:,2)+4*gapy; cnt_cond80(:,4)=80;
% cond 81-96
cnt_cond81=cnt_cond1; cnt_cond81(:,2)=cnt_cond81(:,2)+5*gapy; cnt_cond81(:,4)=81;
cnt_cond82=cnt_cond2; cnt_cond82(:,2)=cnt_cond82(:,2)+5*gapy; cnt_cond82(:,4)=82;
cnt_cond83=cnt_cond3; cnt_cond83(:,2)=cnt_cond83(:,2)+5*gapy; cnt_cond83(:,4)=83;
cnt_cond84=cnt_cond4; cnt_cond84(:,2)=cnt_cond84(:,2)+5*gapy; cnt_cond84(:,4)=84;
cnt_cond85=cnt_cond5; cnt_cond85(:,2)=cnt_cond85(:,2)+5*gapy; cnt_cond85(:,4)=85;
cnt_cond86=cnt_cond6; cnt_cond86(:,2)=cnt_cond86(:,2)+5*gapy; cnt_cond86(:,4)=86;
cnt_cond87=cnt_cond7; cnt_cond87(:,2)=cnt_cond87(:,2)+5*gapy; cnt_cond87(:,4)=87;
cnt_cond88=cnt_cond8; cnt_cond88(:,2)=cnt_cond88(:,2)+5*gapy; cnt_cond88(:,4)=88;
cnt_cond89=cnt_cond9; cnt_cond89(:,2)=cnt_cond89(:,2)+5*gapy; cnt_cond89(:,4)=89;
cnt_cond90=cnt_cond10; cnt_cond90(:,2)=cnt_cond90(:,2)+5*gapy; cnt_cond90(:,4)=90;
cnt_cond91=cnt_cond11; cnt_cond91(:,2)=cnt_cond91(:,2)+5*gapy; cnt_cond91(:,4)=91;
cnt_cond92=cnt_cond12; cnt_cond92(:,2)=cnt_cond92(:,2)+5*gapy; cnt_cond92(:,4)=92;
cnt_cond93=cnt_cond13; cnt_cond93(:,2)=cnt_cond93(:,2)+5*gapy; cnt_cond93(:,4)=93;
cnt_cond94=cnt_cond14; cnt_cond94(:,2)=cnt_cond94(:,2)+5*gapy; cnt_cond94(:,4)=94;
cnt_cond95=cnt_cond15; cnt_cond95(:,2)=cnt_cond95(:,2)+5*gapy; cnt_cond95(:,4)=95;
cnt_cond96=cnt_cond16; cnt_cond96(:,2)=cnt_cond96(:,2)+5*gapy; cnt_cond96(:,4)=96;
% cond 97-112
cnt_cond97=cnt_cond1; cnt_cond97(:,2)=cnt_cond97(:,2)+6*gapy; cnt_cond97(:,4)=97;
cnt_cond98=cnt_cond2; cnt_cond98(:,2)=cnt_cond98(:,2)+6*gapy; cnt_cond98(:,4)=98;
cnt_cond99=cnt_cond3; cnt_cond99(:,2)=cnt_cond99(:,2)+6*gapy; cnt_cond99(:,4)=99;
cnt_cond100=cnt_cond4; cnt_cond100(:,2)=cnt_cond100(:,2)+6*gapy; cnt_cond100(:,4)=100;
cnt_cond101=cnt_cond5; cnt_cond101(:,2)=cnt_cond101(:,2)+6*gapy; cnt_cond101(:,4)=101;
cnt_cond102=cnt_cond6; cnt_cond102(:,2)=cnt_cond102(:,2)+6*gapy; cnt_cond102(:,4)=102;
cnt_cond103=cnt_cond7; cnt_cond103(:,2)=cnt_cond103(:,2)+6*gapy; cnt_cond103(:,4)=103;
cnt_cond104=cnt_cond8; cnt_cond104(:,2)=cnt_cond104(:,2)+6*gapy; cnt_cond104(:,4)=104;
cnt_cond105=cnt_cond9; cnt_cond105(:,2)=cnt_cond105(:,2)+6*gapy; cnt_cond105(:,4)=105;
cnt_cond106=cnt_cond10; cnt_cond106(:,2)=cnt_cond106(:,2)+6*gapy; cnt_cond106(:,4)=106;
cnt_cond107=cnt_cond11; cnt_cond107(:,2)=cnt_cond107(:,2)+6*gapy; cnt_cond107(:,4)=107;
cnt_cond108=cnt_cond12; cnt_cond108(:,2)=cnt_cond108(:,2)+6*gapy; cnt_cond108(:,4)=108;
cnt_cond109=cnt_cond13; cnt_cond109(:,2)=cnt_cond109(:,2)+6*gapy; cnt_cond109(:,4)=109;
cnt_cond110=cnt_cond14; cnt_cond110(:,2)=cnt_cond110(:,2)+6*gapy; cnt_cond110(:,4)=110;
cnt_cond111=cnt_cond15; cnt_cond111(:,2)=cnt_cond111(:,2)+6*gapy; cnt_cond111(:,4)=111;
cnt_cond112=cnt_cond16; cnt_cond112(:,2)=cnt_cond112(:,2)+6*gapy; cnt_cond112(:,4)=112;
% cond 113-128
cnt_cond113=cnt_cond1; cnt_cond113(:,2)=cnt_cond113(:,2)+7*gapy; cnt_cond113(:,4)=113;
cnt_cond114=cnt_cond2; cnt_cond114(:,2)=cnt_cond114(:,2)+7*gapy; cnt_cond114(:,4)=114;
cnt_cond115=cnt_cond3; cnt_cond115(:,2)=cnt_cond115(:,2)+7*gapy; cnt_cond115(:,4)=115;
cnt_cond116=cnt_cond4; cnt_cond116(:,2)=cnt_cond116(:,2)+7*gapy; cnt_cond116(:,4)=116;
cnt_cond117=cnt_cond5; cnt_cond117(:,2)=cnt_cond117(:,2)+7*gapy; cnt_cond117(:,4)=117;
cnt_cond118=cnt_cond6; cnt_cond118(:,2)=cnt_cond118(:,2)+7*gapy; cnt_cond118(:,4)=118;
cnt_cond119=cnt_cond7; cnt_cond119(:,2)=cnt_cond119(:,2)+7*gapy; cnt_cond119(:,4)=119;
cnt_cond120=cnt_cond8; cnt_cond120(:,2)=cnt_cond120(:,2)+7*gapy; cnt_cond120(:,4)=120;
cnt_cond121=cnt_cond9; cnt_cond121(:,2)=cnt_cond121(:,2)+7*gapy; cnt_cond121(:,4)=121;
cnt_cond122=cnt_cond10; cnt_cond122(:,2)=cnt_cond122(:,2)+7*gapy; cnt_cond122(:,4)=122;
cnt_cond123=cnt_cond11; cnt_cond123(:,2)=cnt_cond123(:,2)+7*gapy; cnt_cond123(:,4)=123;
cnt_cond124=cnt_cond12; cnt_cond124(:,2)=cnt_cond124(:,2)+7*gapy; cnt_cond124(:,4)=124;
cnt_cond125=cnt_cond13; cnt_cond125(:,2)=cnt_cond125(:,2)+7*gapy; cnt_cond125(:,4)=125;
cnt_cond126=cnt_cond14; cnt_cond126(:,2)=cnt_cond126(:,2)+7*gapy; cnt_cond126(:,4)=126;
cnt_cond127=cnt_cond15; cnt_cond127(:,2)=cnt_cond127(:,2)+7*gapy; cnt_cond127(:,4)=127;
cnt_cond128=cnt_cond16; cnt_cond128(:,2)=cnt_cond128(:,2)+7*gapy; cnt_cond128(:,4)=128;


Cnt = [cnt_cond1;cnt_cond2;cnt_cond3;cnt_cond4;cnt_cond5;cnt_cond6;cnt_cond7;cnt_cond8;cnt_cond9;cnt_cond10;cnt_cond11;cnt_cond12; ...
    cnt_cond13;cnt_cond14;cnt_cond15;cnt_cond16;cnt_cond17;cnt_cond18;cnt_cond19;cnt_cond20;cnt_cond21;cnt_cond22;cnt_cond23;cnt_cond24; ...
    cnt_cond25;cnt_cond26;cnt_cond27;cnt_cond28;cnt_cond29;cnt_cond30;cnt_cond31;cnt_cond32;cnt_cond33;cnt_cond34;cnt_cond35;cnt_cond36; ...
    cnt_cond37;cnt_cond38;cnt_cond39;cnt_cond40;cnt_cond41;cnt_cond42;cnt_cond43;cnt_cond44;cnt_cond45;cnt_cond46;cnt_cond47;cnt_cond48; ...
    cnt_cond49;cnt_cond50;cnt_cond51;cnt_cond52;cnt_cond53;cnt_cond54;cnt_cond55;cnt_cond56;cnt_cond57;cnt_cond58;cnt_cond59;cnt_cond60; ...
    cnt_cond61;cnt_cond62;cnt_cond63;cnt_cond64;cnt_cond65;cnt_cond66;cnt_cond67;cnt_cond68;cnt_cond69;cnt_cond70;cnt_cond71;cnt_cond72; ...
    cnt_cond73;cnt_cond74;cnt_cond75;cnt_cond76;cnt_cond77;cnt_cond78;cnt_cond79;cnt_cond80;cnt_cond81;cnt_cond82;cnt_cond83;cnt_cond84; ...
    cnt_cond85;cnt_cond86;cnt_cond87;cnt_cond88;cnt_cond89;cnt_cond90;cnt_cond91;cnt_cond92;cnt_cond93;cnt_cond94;cnt_cond95;cnt_cond96; ...
    cnt_cond97;cnt_cond98;cnt_cond99;cnt_cond100;cnt_cond101;cnt_cond102;cnt_cond103;cnt_cond104;cnt_cond105;cnt_cond106;cnt_cond107;cnt_cond108; ...
    cnt_cond109;cnt_cond110;cnt_cond111;cnt_cond112;cnt_cond113;cnt_cond114;cnt_cond115;cnt_cond116;cnt_cond117;cnt_cond118;cnt_cond119; ...
    cnt_cond120;cnt_cond121;cnt_cond122;cnt_cond123;cnt_cond124;cnt_cond125;cnt_cond126;cnt_cond127;cnt_cond128]; % centers of conductors
clear cnt_cond1 cnt_cond2 cnt_cond3 cnt_cond4 cnt_cond5 cnt_cond6 cnt_cond7 cnt_cond8 cnt_cond9 cnt_cond10 cnt_cond11 cnt_cond12 ...
    cnt_cond13 cnt_cond14 cnt_cond15 cnt_cond16 cnt_cond17 cnt_cond18 cnt_cond19 cnt_cond20 cnt_cond21 cnt_cond22 cnt_cond23 cnt_cond24 ...
    cnt_cond25 cnt_cond26 cnt_cond27 cnt_cond28 cnt_cond29 cnt_cond30 cnt_cond31 cnt_cond32 cnt_cond33 cnt_cond34 cnt_cond35 cnt_cond36  ...
    cnt_cond37 cnt_cond38 cnt_cond39 cnt_cond40 cnt_cond41 cnt_cond42 cnt_cond43 cnt_cond44 cnt_cond45 cnt_cond46 cnt_cond47 cnt_cond48  ...
    cnt_cond49 cnt_cond50 cnt_cond51 cnt_cond52 cnt_cond53 cnt_cond54 cnt_cond55 cnt_cond56 cnt_cond57 cnt_cond58 cnt_cond59 cnt_cond60  ...
    cnt_cond61 cnt_cond62 cnt_cond63 cnt_cond64 cnt_cond65 cnt_cond66 cnt_cond67 cnt_cond68 cnt_cond69 cnt_cond70 cnt_cond71 cnt_cond72  ...
    cnt_cond73 cnt_cond74 cnt_cond75 cnt_cond76 cnt_cond77 cnt_cond78 cnt_cond79 cnt_cond80 cnt_cond81 cnt_cond82 cnt_cond83 cnt_cond84  ...
    cnt_cond85 cnt_cond86 cnt_cond87 cnt_cond88 cnt_cond89 cnt_cond90 cnt_cond91 cnt_cond92 cnt_cond93 cnt_cond94 cnt_cond95 cnt_cond96  ...
    cnt_cond97 cnt_cond98cnt_cond99 cnt_cond100 cnt_cond101 cnt_cond102 cnt_cond103 cnt_cond104 cnt_cond105 cnt_cond106 cnt_cond107 cnt_cond108  ...
    cnt_cond109 cnt_cond110 cnt_cond111 cnt_cond112 cnt_cond113 cnt_cond114 cnt_cond115 cnt_cond16 cnt_cond117 cnt_cond118 cnt_cond119  ...
    cnt_cond120 cnt_cond121 cnt_cond122 cnt_cond123 cnt_cond124 cnt_cond125 cnt_cond126 cnt_cond127 cnt_cond128

Dim1=[len1 wid hei]; Dim2=[len2 wid hei]; Dim3=[len3 wid hei]; Dim4=[len4 wid hei]; Dim5=[len5 wid hei]; Dim6=[len6 wid hei];
Dim7=[len7 wid hei]; Dim8=[len8 wid hei]; Dim9=[len9 wid hei]; Dim10=[len10 wid hei]; Dim11=[len11 wid hei]; Dim12=[len12 wid hei];
Dim_temp=[Dim1;Dim2;Dim3;Dim4;Dim5;Dim6;Dim7;Dim8;Dim9;Dim10;Dim11;Dim12];
Dims=[];
for ii=1:128
    Dims=[Dims; Dim_temp];
end

Orients_tmp=['x';'y';'x';'y';'x';'y';'x';'y';'x';'y';'x';'y']; % orientations of conductors %
Orients=[];
for ii=1:128
    Orients=[Orients; Orients_tmp];
end
clear Dim1 Dim2 Dim3 Dim4 Dim5 Dim6 Dim7 Dim8 Dim9 Dim10 Dim11 Dim12 Orients_tmp Dim_temp


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
bbox_max=[gapx*15+len1 gapy*7+len2+wid hei]+1e-12; % max coordinates of bbox

% -------------------------------------------------------------------------
%                  Input for Ports
% -------------------------------------------------------------------------
% At the end of this part, we need structures pnt_lft{xx} and pnt_rght{xx} which
% contains the coordinates of nodes on both sides of xxth port
port_direction=[1 2];%normal direction of left and right port. 1-x 2-y 3-z
% defining the nodes in first port
pnt_lft=cell(num_ports,1); %exc
pnt_rght=cell(num_ports,1); %gnd
pnt_lft{1}=zeros(round(wid/Res)*round(hei/Res),4);
pnt_rght{1}=zeros(round(wid/Res)*round(hei/Res),4);
dum=1;

for kk=1:round(wid/Res)
    for ll=1:round(hei/Res)
        pnt_rght{1}(dum,1:4)=[(2*kk-1)*(0.5*Res)+2*wid+3*spa 3*spa+3*wid (2*ll-1)*(0.5*Res) 2];
        pnt_lft{1}(dum,1:4)=[0 (2*kk-1)*(0.5*Res) (2*ll-1)*(0.5*Res) 1];
        dum=dum+1;
    end
end
% defining the nodes in remaining ports
for kk=2:16
    pnt_lft{kk}(:,1)=pnt_lft{1}(:,1) + gapx*(kk-1);
    pnt_lft{kk}(:,[2 3 4])=pnt_lft{1}(:,[2 3 4]);
    pnt_rght{kk}(:,1)=pnt_rght{1}(:,1) + gapx*(kk-1);
    pnt_rght{kk}(:,[2 3 4])=pnt_rght{1}(:,[2 3 4]);
end
for kk=17:32
    pnt_lft{kk}(:,2)=pnt_lft{kk-16}(:,2) + gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-16}(:,2) + gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-16}(:,[1 3 4]);
end
for kk=33:48
    pnt_lft{kk}(:,2)=pnt_lft{kk-2*16}(:,2) + 2*gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-2*16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-2*16}(:,2) + 2*gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-2*16}(:,[1 3 4]);
end
for kk=49:64
    pnt_lft{kk}(:,2)=pnt_lft{kk-3*16}(:,2) + 3*gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-3*16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-3*16}(:,2) + 3*gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-3*16}(:,[1 3 4]);
end
for kk=65:80
    pnt_lft{kk}(:,2)=pnt_lft{kk-4*16}(:,2) + 4*gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-4*16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-4*16}(:,2) + 4*gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-4*16}(:,[1 3 4]);
end
for kk=81:96
    pnt_lft{kk}(:,2)=pnt_lft{kk-5*16}(:,2) + 5*gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-5*16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-5*16}(:,2) + 5*gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-5*16}(:,[1 3 4]);
end
for kk=97:112
    pnt_lft{kk}(:,2)=pnt_lft{kk-6*16}(:,2) + 6*gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-6*16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-6*16}(:,2) + 6*gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-6*16}(:,[1 3 4]);
end
for kk=113:128
    pnt_lft{kk}(:,2)=pnt_lft{kk-7*16}(:,2) + 7*gapy;
    pnt_lft{kk}(:,[1 3 4])=pnt_lft{kk-7*16}(:,[1 3 4]);
    pnt_rght{kk}(:,2)=pnt_rght{kk-7*16}(:,2) + 7*gapy;
    pnt_rght{kk}(:,[1 3 4])=pnt_rght{kk-7*16}(:,[1 3 4]);
end
% defining nodes connected ground if conductors without ports exist; if
% there is no, then leave as a empty array.

%pnt_well_cond=[pnt_lft{1}(1,1) pnt_lft{1}(1,2) + dist_btw_conds pnt_lft{1}(1,3);];
pnt_well_cond=[];

