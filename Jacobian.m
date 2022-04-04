function [J,B,Ro_ab,Rg_ab,Rw_ab]=Jacobian(index,numi,numj,Sw_new,Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,count,PVT_OIL,PVT_WATER,PVT_GAS,OW,OG,CR)
% This function is used to calculate the Jacobian matrix and the residual

J=zeros(3*count);
B=zeros(3*count,1);
rhoo=zeros(numi,numj);
Bo=zeros(numi,numj);
muo=zeros(numi,numj);
Rso=zeros(numi,numj);
krw=zeros(numi,numj);
krg=zeros(numi,numj);
kro=zeros(numi,numj);
Pcow=zeros(numi,numj);
Pcgo=zeros(numi,numj);
rhow=zeros(numi,numj);
Bw=zeros(numi,numj);
muw=zeros(numi,numj);
rhog=zeros(numi,numj);
Bg=zeros(numi,numj);
mug=zeros(numi,numj);
for i=1:numi
    for j=1:numj        
        if index(i,j)~=0
            [rhoo(i,j),Bo(i,j),muo(i,j),Rso(i,j)]=PVT_oil(Po(i,j),PVT_OIL);
            [krw(i,j),krg(i,j),kro(i,j),Pcow(i,j),Pcgo(i,j)]=relaperm(Sw_new(i,j),Sg_new(i,j),OW,OG);
            [rhow(i,j),Bw(i,j),muw(i,j)]=PVT_water(Po(i,j)-Pcow(i,j),PVT_WATER);
            [rhog(i,j),Bg(i,j),mug(i,j)]=PVT_gas(Po(i,j)+Pcgo(i,j),PVT_GAS);
        end
    end
end   
for num=1:count
    [a,b]=find(grid_num==num);
% I. LHS, the Jacobian Matrix   
% 1. Self    
    num_target=num;
        [Jcell,Ro_ab,Rw_ab,Rg_ab]=Jacobian_cell(num_target,index,a,b,numi,numj,...
            Sw_new,Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,...
            Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,PVT_OIL,PVT_WATER,PVT_GAS,OW,OG,CR...
            ,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
    for i=-2:0
        for j=-2:0
            J(3*num+i,3*num_target+j)=Jcell(i+3,j+3);
        end
    end

% 2. North

    if index(a+1,b)~=0 % North is not boundary
        num_target=grid_num(a+1,b);
        [Jcell,~,~,~]=Jacobian_cell(num_target,index,a,b,numi,numj,Sw_new,Sg_new,...
            Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,...
            phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,PVT_OIL,PVT_WATER,PVT_GAS,OW,OG,CR...
            ,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
    
%     for i=-2:0
%         for j=-2:0
%             J(3*(num+1)+i,3*num+j)=Jcell(i+3,j+3);
%         end
%     end
     for i=-2:0
        for j=-2:0
            J(3*(num)+i,3*num_target+j)=Jcell(i+3,j+3);
        end
    end
    end
    
% 3. South

    if index(a-1,b)~=0 % South is not boundary
        num_target=grid_num(a-1,b);
        [Jcell,~,~,~]=Jacobian_cell(num_target,index,a,b,numi,numj,Sw_new,...
            Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,...
            Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,PVT_OIL,...
            PVT_WATER,PVT_GAS,OW,OG,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
   
%     for i=-2:0
%         for j=-2:0
%             J(3*(num-1)+i,3*num+j)=Jcell(i+3,j+3);
%         end
%     end
    for i=-2:0
        for j=-2:0
            J(3*(num)+i,3*num_target+j)=Jcell(i+3,j+3);
        end
    end
    end
% 4. East

    if index(a,b+1)~=0 % East is not boundary
        num_target=grid_num(a,b+1);
        [Jcell,~,~,~]=Jacobian_cell(num_target,index,a,b,numi,numj,Sw_new,...
            Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,...
            Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,PVT_OIL,PVT_WATER,...
            PVT_GAS,OW,OG,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
  
    for i=-2:0
        for j=-2:0
            J(3*num+i,3*num_target+j)=Jcell(i+3,j+3);
        end
    end
 
    end

% 5. West

    if index(a,b-1)~=0 % West is not boundary
        num_target=grid_num(a,b-1);
        [Jcell,~,~,~]=Jacobian_cell(num_target,index,a,b,numi,numj,Sw_new,...
            Sg_new,Sw_old,Sg_old,Po,G,dx,dy,dz,phi_old,Bo_old,Bg_old,Bw_old,...
            Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,grid_num,dt,kx,ky,PVT_OIL,PVT_WATER,...
            PVT_GAS,OW,OG,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso);
    
    for i=-2:0
        for j=-2:0
            J(3*num+i,3*num_target+j)=Jcell(i+3,j+3);
        end
    end
    end

% II. RHS, the residual term
    
   B(3*num-2)=-Ro_ab;
   B(3*num-1)=-Rw_ab;
   B(3*num)=-Rg_ab;
end
end