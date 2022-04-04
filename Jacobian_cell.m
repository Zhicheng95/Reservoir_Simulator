function [Jcell,Ro_ab,Rw_ab,Rg_ab]=Jacobian_cell(num_target,index,a,b,numi,numj,...
    Sw_new,Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,Rso_old,...
    Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,PVT_OIL,PVT_WATER,PVT_GAS,OW,OG,CR...
    ,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso)
% Derivative of grid(a,b) with grid block(number) given a stir
ep=1e-6;% epslon of pressrue
es=0.0000001;% epslon of saturation
Jcell=zeros(3);
% rhoo=zeros(numi,numj);
% Bo=zeros(numi,numj);
% muo=zeros(numi,numj);
% Rso=zeros(numi,numj);
% krw=zeros(numi,numj);
% krg=zeros(numi,numj);
% kro=zeros(numi,numj);
% Pcow=zeros(numi,numj);
% Pcgo=zeros(numi,numj);
% rhow=zeros(numi,numj);
% Bw=zeros(numi,numj);
% muw=zeros(numi,numj);
% rhog=zeros(numi,numj);
% Bg=zeros(numi,numj);
% mug=zeros(numi,numj);
% for i=1:numi
%     for j=1:numj        
%         if index(i,j)~=0
%             [rhoo(i,j),Bo(i,j),muo(i,j),Rso(i,j)]=PVT_oil(Po(i,j),PVT_OIL);
%             [krw(i,j),krg(i,j),kro(i,j),Pcow(i,j),Pcgo(i,j)]=relaperm(Sw_new(i,j),Sg_new(i,j),OW,OG);
%             [rhow(i,j),Bw(i,j),muw(i,j)]=PVT_water(Po(i,j)-Pcow(i,j),PVT_WATER);
%             [rhog(i,j),Bg(i,j),mug(i,j)]=PVT_gas(Po(i,j)+Pcgo(i,j),PVT_GAS);
%         end
%     end
% end       
[Ro_ini,Rw_ini,Rg_ini]=Residual3(Sw_new,Sg_new,Sw_old,Sg_old,Po,numi,numj,index,G,dx,dy,dz,phi_old,...
    Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,kx,ky,dt,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
for i=1:numi
    for j=1:numj
        if grid_num(i,j)==num_target
            x=i;y=j;
        end
    end
end

Ro_ab=Ro_ini(a,b);
Rw_ab=Rw_ini(a,b);
Rg_ab=Rg_ini(a,b);
%% Derivative of Ro,Rw and Rg respect to the grid(x,y) with Pressure change,
% of grid(a,b)
Po(x,y)=Po(x,y)+ep;

[rhoo(x,y),Bo(x,y),muo(x,y),Rso(x,y)]=PVT_oil(Po(x,y),PVT_OIL);
[krw(x,y),krg(x,y),kro(x,y),Pcow(x,y),Pcgo(x,y)]=relaperm(Sw_new(x,y),Sg_new(x,y),OW,OG);
[rhow(x,y),Bw(x,y),muw(x,y)]=PVT_water(Po(x,y)-Pcow(x,y),PVT_WATER);
[rhog(x,y),Bg(x,y),mug(x,y)]=PVT_gas(Po(x,y)+Pcgo(x,y),PVT_GAS);

[Ro_p,Rw_p,Rg_p]=Residual3(Sw_new,Sg_new,Sw_old,Sg_old,Po,numi,numj,index,G,dx,dy,dz,phi_old,...
    Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,kx,ky,dt,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
Jcell(1,1)=(Ro_p(a,b)-Ro_ini(a,b))/ep;
Jcell(2,1)=(Rw_p(a,b)-Rw_ini(a,b))/ep;
Jcell(3,1)=(Rg_p(a,b)-Rg_ini(a,b))/ep;
Po(x,y)=Po(x,y)-ep;

[rhoo(x,y),Bo(x,y),muo(x,y),Rso(x,y)]=PVT_oil(Po(x,y),PVT_OIL);
[krw(x,y),krg(x,y),kro(x,y),Pcow(x,y),Pcgo(x,y)]=relaperm(Sw_new(x,y),Sg_new(x,y),OW,OG);
[rhow(x,y),Bw(x,y),muw(x,y)]=PVT_water(Po(x,y)-Pcow(x,y),PVT_WATER);
[rhog(x,y),Bg(x,y),mug(x,y)]=PVT_gas(Po(x,y)+Pcgo(x,y),PVT_GAS);

% Derivative of Ro,Rw and Rg respect to the grid(x,y) with Sw change,
% of grid(a,b)
Sw_new(x,y)=Sw_new(x,y)+es;

[rhoo(x,y),Bo(x,y),muo(x,y),Rso(x,y)]=PVT_oil(Po(x,y),PVT_OIL);
[krw(x,y),krg(x,y),kro(x,y),Pcow(x,y),Pcgo(x,y)]=relaperm(Sw_new(x,y),Sg_new(x,y),OW,OG);
[rhow(x,y),Bw(x,y),muw(x,y)]=PVT_water(Po(x,y)-Pcow(x,y),PVT_WATER);
[rhog(x,y),Bg(x,y),mug(x,y)]=PVT_gas(Po(x,y)+Pcgo(x,y),PVT_GAS);

[Ro_w,Rw_w,Rg_w]=Residual3(Sw_new,Sg_new,Sw_old,Sg_old,Po,numi,numj,index,G,dx,dy,dz,phi_old,...
    Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,kx,ky,dt,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
Jcell(1,2)=(Ro_w(a,b)-Ro_ini(a,b))/es;
Jcell(2,2)=(Rw_w(a,b)-Rw_ini(a,b))/es;
Jcell(3,2)=(Rg_w(a,b)-Rg_ini(a,b))/es;
Sw_new(x,y)=Sw_new(x,y)-es;

[rhoo(x,y),Bo(x,y),muo(x,y),Rso(x,y)]=PVT_oil(Po(x,y),PVT_OIL);
[krw(x,y),krg(x,y),kro(x,y),Pcow(x,y),Pcgo(x,y)]=relaperm(Sw_new(x,y),Sg_new(x,y),OW,OG);
[rhow(x,y),Bw(x,y),muw(x,y)]=PVT_water(Po(x,y)-Pcow(x,y),PVT_WATER);
[rhog(x,y),Bg(x,y),mug(x,y)]=PVT_gas(Po(x,y)+Pcgo(x,y),PVT_GAS);

% Derivative of Ro,Rw and Rg respect to the grid(x,y) with Sg change,
% of grid(a,b)
Sg_new(x,y)=Sg_new(x,y)+es;

[rhoo(x,y),Bo(x,y),muo(x,y),Rso(x,y)]=PVT_oil(Po(x,y),PVT_OIL);
[krw(x,y),krg(x,y),kro(x,y),Pcow(x,y),Pcgo(x,y)]=relaperm(Sw_new(x,y),Sg_new(x,y),OW,OG);
[rhow(x,y),Bw(x,y),muw(x,y)]=PVT_water(Po(x,y)-Pcow(x,y),PVT_WATER);
[rhog(x,y),Bg(x,y),mug(x,y)]=PVT_gas(Po(x,y)+Pcgo(x,y),PVT_GAS);

[Ro_g,Rw_g,Rg_g]=Residual3(Sw_new,Sg_new,Sw_old,Sg_old,Po,numi,numj,index,G,dx,dy,dz,phi_old,...
    Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,kx,ky,dt,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
Jcell(1,3)=(Ro_g(a,b)-Ro_ini(a,b))/es;
Jcell(2,3)=(Rw_g(a,b)-Rw_ini(a,b))/es;
Jcell(3,3)=(Rg_g(a,b)-Rg_ini(a,b))/es;
Sg_new(x,y)=Sg_new(x,y)-es;
% [rhoo(x,y),Bo(x,y),muo(x,y),Rso(x,y)]=PVT_oil(Po(x,y),PVT_OIL);
% [krw(x,y),krg(x,y),kro(x,y),Pcow(x,y),Pcgo(x,y)]=relaperm(Sw(x,y),Sg(x,y),OW,OG);
% [rhow(x,y),Bw(x,y),muw(x,y)]=PVT_water(Po(x,y)-Pcow(x,y),PVT_WATER);
% [rhog(x,y),Bg(x,y),mug(x,y)]=PVT_gas(Po(x,y)+Pcgo(x,y),PVT_GAS);
end
