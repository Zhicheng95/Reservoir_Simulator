clc
clear

%% This is the main function for the sample Jacobian and calculation for the sample case

index=load('index.txt');

dx=load('DX.txt');
dy=load('DY.txt');
dz=load('DZ.txt');
phi_ini=load('Porosity.txt');
kx=load('perm_NS.txt');
G=load('G.txt');

[index,grid_num]=matrix_tp(index,1);
[dx,~]=matrix_tp(dx,0);
[dy,~]=matrix_tp(dy,0);
[dz,~]=matrix_tp(dz,0);
[phi_ini,~]=matrix_tp(phi_ini,0);
[kx,~]=matrix_tp(kx,0);
ky=0.5*kx;
[G,~]=matrix_tp(G,0);


[numi,numj]=size(index);

CR=3*10^(-6); % psi^-1


count=1;
for j=1:numj
    for i=1:numi
        if index(i,j)~=0
            grid_num(i,j)=count;count=count+1;
        end
        
    end
end
count=count-1;

Sw_ini=zeros(numi,numj);
Sg_ini=zeros(numi,numj);
Po_ini=zeros(numi,numj);

for i=1:numi
    for j=1:numj
        if index(i,j)~=0
            Sw_ini(i,j)=0.3;
            Sg_ini(i,j)=0.1;
            Po_ini(i,j)=2000;
        end
    end
end

dt=5;
PVT_OIL=load('PVT_oil.txt');
PVT_GAS=load('PVT_gas.txt');
PVT_WATER=load('PVT_water.txt');
OW=load('OW relative perm.txt');
OG=load('OG relative perm.txt');

sp_condition=zeros(numi,numj);
sp_condition(46,22)=8;
sp_condition(42,16)=2;
sp_condition(40,10)=8;
sp_condition(37,23)=8;
sp_condition(34,12)=1;
sp_condition(31,15)=5;
sp_condition(31,12)=8;
sp_condition(30,6)=4;
sp_condition(26,19)=8;
sp_condition(26,7)=8;

sp_data=zeros(numi,numj);
sp_data(46,22)=1600;
sp_data(42,16)=-25;
sp_data(40,10)=1600;
sp_data(37,23)=1600;
sp_data(34,12)=-50000;
sp_data(31,15)=-30;
sp_data(31,12)=1600;
sp_data(30,6)=-900000;
sp_data(26,19)=1600;
sp_data(26,7)=1600;

% sp_condition=zeros(numi,numj);
% sp_condition(46,22)=2;
% sp_condition(42,16)=2;
% sp_condition(40,10)=2;
% sp_condition(37,23)=2;
% sp_condition(34,12)=2;
% sp_condition(31,15)=2;
% sp_condition(31,12)=2;
% sp_condition(30,6)=2;
% sp_condition(26,19)=2;
% sp_condition(26,7)=2;
% 
% sp_data=zeros(numi,numj);
% sp_data(46,22)=-25;
% sp_data(42,16)=-25;
% sp_data(40,10)=-25;
% sp_data(37,23)=-25;
% sp_data(34,12)=-25;
% sp_data(31,15)=-25;
% sp_data(31,12)=-25;
% sp_data(30,6)=-25;
% sp_data(26,19)=-25;
% sp_data(26,7)=-25;

Sw_new=Sw_ini;
Sg_new=Sg_ini;
Sw_old=Sw_ini;
Sg_old=Sg_ini;
Po=Po_ini;
Po_old=Po;
phi_old=phi_ini;

tag_o=0;
tag_w=0;
dPo=zeros(count,1);
dSw=zeros(count,1);
dSg=zeros(count,1);


Bo_old=zeros(numi,numj);
Rso_old=zeros(numi,numj);
Bw_old=zeros(numi,numj);
Bg_old=zeros(numi,numj);
Pcow_old=zeros(numi,numj);
Pcgo_old=zeros(numi,numj);

% For first iteration
for i=1:numi
    for j=1:numj
        if index(i,j)~=0
            [~,~,~,Pcow_old(i,j),Pcgo_old(i,j)]=relaperm(Sw_old(i,j),Sg_old(i,j),OW,OG);
            [~,Bo_old(i,j),~,Rso_old(i,j)]=PVT_oil(Po_old(i,j),PVT_OIL);
            [~,Bw_old(i,j),~]=PVT_water(Po_old(i,j)-Pcow_old(i,j),PVT_WATER);
            [~,Bg_old(i,j),~]=PVT_gas(Po_old(i,j)+Pcgo_old(i,j),PVT_GAS);
        end
    end
end
Bo_ini=Bo_old;
Bw_ini=Bw_old;
Bg_ini=Bg_old;

numt=360/dt;
Pov=zeros(numi,numj,numt+1);
Pov(:,:,1)=Po;
Qov=zeros(numi,numj,numt+1);
Qwv=zeros(numi,numj,numt+1);
Qgv=zeros(numi,numj,numt+1);
Swv=zeros(numi,numj,numt+1);
Sgv=zeros(numi,numj,numt+1);

[Qo,Qw,Qg,~,~]=flowrate(sp_condition,sp_data,PVT_WATER,PVT_OIL,PVT_GAS,OW,OG,numi,numj,Po,Sw_new,Sg_new,index,kx,ky,dx,dy,dz);
Qov(:,:,1)=Qo;    
Qwv(:,:,1)=Qw;  
Qgv(:,:,1)=Qg;

Sw_old=Sw_new;
Sg_old=Sg_new;
Po_old=Po;
iterv=zeros(numt+1,1);

oil_inc=zeros(numt,1);
gas_inc=zeros(numt,1);
water_inc=zeros(numt,1);
oil_com=zeros(numt,1);
gas_com=zeros(numt,1);
water_com=zeros(numt,1);
Ro_ab=zeros(numt,1);
Rg_ab=zeros(numt,1);
Rw_ab=zeros(numt,1);

%%
for time=2:numt+1
iter=1;
tag_inc=0;
tag_com=0;
while ((tag_inc==0 || tag_com==0) && iter<=10) || max(B)>0.01 || iter==1
[Qo,Qw,Qg,~,~]=flowrate(sp_condition,sp_data,PVT_WATER,PVT_OIL,PVT_GAS,OW,OG,numi,numj,Po,Sw_new,Sg_new,index,kx,ky,dx,dy,dz);  
[J,B,Ro_ab(time-1),Rg_ab(time-1),Rw_ab(time-1)]=Jacobian(index,numi,numj,Sw_new,Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,count,PVT_OIL,PVT_WATER,PVT_GAS,OW,OG,CR);

X=J\B;

for i=1:count
    dPo(i)=X(3*i-2);
    dSw(i)=X(3*i-1);
    dSg(i)=X(3*i);
    [a,b]=find(grid_num==i);
    Po(a,b)=Po(a,b)+dPo(i);
    Sw_new(a,b)=Sw_new(a,b)+dSw(i);
    Sg_new(a,b)=Sg_new(a,b)+dSg(i);
end
[Qo,Qw,Qg,~,~]=flowrate(sp_condition,sp_data,PVT_WATER,PVT_OIL,PVT_GAS,OW,OG,numi,numj,Po,Sw_new,Sg_new,index,kx,ky,dx,dy,dz);
Qov(:,:,time)=Qo;       
Qwv(:,:,time)=Qw;
Qgv(:,:,time)=Qg;
[tag_inc,tag_com]=MBcheck(numi,numj,numt,time,Po,Po_old,Sw_old,Sg_old,Po_ini,phi_ini,...
    Sw_ini,Sg_ini,Sw_new,Sg_new,Qo,Qw,Qg,Qov,Qwv,Qgv,PVT_GAS,PVT_WATER,PVT_OIL,OW,OG,CR,dt,index,dx,dy,dz,Bo_ini,Bw_ini,Bg_ini);

iter=iter+1;

end
[~,~,oil_inc(time),water_inc(time),gas_inc(time),oil_com(time),water_com(time),gas_com(time)]=MBcheck(numi,numj,numt,time,Po,Po_old,Sw_old,Sg_old,Po_ini,phi_ini,...
    Sw_ini,Sg_ini,Sw_new,Sg_new,Qo,Qw,Qg,Qov,Qwv,Qgv,PVT_GAS,PVT_WATER,PVT_OIL,OW,OG,CR,dt,index,dx,dy,dz,Bo_ini,Bw_ini,Bg_ini);
Pov(:,:,time)=Po;
Qov(:,:,time)=Qo;
Qwv(:,:,time)=Qw;
Qgv(:,:,time)=Qg;
Swv(:,:,time)=Sw_new;
Sgv(:,:,time)=Sg_new;
iterv(time)=iter;
Sw_old=Sw_new;
Sg_old=Sg_new;
Po_old=Po;
for i=1:numi 
    for j=1:numj
        if index(i,j)~=0
            [~,~,~,Pcow_old(i,j),Pcgo_old(i,j)]=relaperm(Sw_old(i,j),Sg_old(i,j),OW,OG);
            phi_old(i,j)=phi_ini(i,j)*(1+CR*(Po_old(i,j)-Po_ini(i,j)));
            [~,Bo_old(i,j),~,Rso_old(i,j)]=PVT_oil(Po_old(i,j),PVT_OIL);
            [~,Bw_old(i,j),~]=PVT_water(Po_old(i,j)-Pcow_old(i,j),PVT_WATER);
            [~,Bg_old(i,j),~]=PVT_gas(Po_old(i,j)+Pcgo_old(i,j),PVT_GAS);
        end
    end
end

end

