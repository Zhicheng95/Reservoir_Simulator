function [Ro,Rw,Rg]=Residual3(Sw_new,Sg_new,Sw_old,Sg_old,Po,numi,numj,index,G,dx,dy,dz,phi_old,...
    Bo_old,Bg_old,Bw_old,Rso_old,Po_ini,phi_ini,Qo,Qw,Qg,kx,ky,dt,CR,rhoo,Bo,muo,krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug,Rso)
% This function is used to calcualte the residual for the 3 phases systems
% for the 2D reservoir
% These transmisibility terms are based on new Sw,Sg iteration value

% Sw_new, Sg_new is used to note down the n+1 time step value;
% Sw_old, Sg_old is used to note down the n time step value;

% CR=0*10^(-6);
% rhoo=zeros(numi,numj);
% Bo=zeros(numi,numj);
% rhow=zeros(numi,numj);
% Bw=zeros(numi,numj);
% rhog=zeros(numi,numj);
% Bg=zeros(numi,numj);
% Pcow=zeros(numi,numj);
% Pcgo=zeros(numi,numj);
phi_new=zeros(numi,numj);
% Rso=zeros(numi,numj);

Ro=zeros(numi,numj);
Rw=zeros(numi,numj);
Rg=zeros(numi,numj);
Pw=zeros(numi,numj);
Pg=zeros(numi,numj);



for i=1:numi
    for j=1:numj
        if index(i,j)~=0
%             [rhoo(i,j),Bo(i,j),~,Rso(i,j)]=PVT_oil(Po(i,j),PVT_OIL);
%             [~,~,~,Pcow(i,j),Pcgo(i,j)]=relaperm(Sw_new(i,j),Sg_new(i,j),OW,OG);
%             [rhow(i,j),Bw(i,j),~]=PVT_water(Po(i,j)-Pcow(i,j),PVT_WATER);
%             [rhog(i,j),Bg(i,j),~]=PVT_gas(Po(i,j)+Pcgo(i,j),PVT_GAS);
%             
            phi_new(i,j)=phi_ini(i,j)*(1+CR*(Po(i,j)-Po_ini(i,j)));
            Pw(i,j)=Po(i,j)-Pcow(i,j);
            Pg(i,j)=Po(i,j)+Pcgo(i,j);
        end
    end
end

[TN_o,TN_w,TN_g,TS_o,TS_w,TS_g,TE_o,TE_w,TE_g,TW_o,TW_w,TW_g]=Transmisibility3(Po,numi,numj,index,G,dx,dy,dz,kx,ky,rhoo,Bo,muo,...
   krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug);

for i=1:numi
    for j=1:numj
        if index(i,j)~=0
            Ro(i,j)=TN_o(i,j)*(Po(i+1,j)-Po(i,j)-1/144*(1/2*(rhoo(i,j)+rhoo(i+1,j)))*(G(i+1,j)-G(i,j)))-...
                    TS_o(i,j)*(Po(i,j)-Po(i-1,j)-1/144*(1/2*(rhoo(i,j)+rhoo(i-1,j)))*(G(i,j)-G(i-1,j)))+...
                    TE_o(i,j)*(Po(i,j+1)-Po(i,j)-1/144*(1/2*(rhoo(i,j)+rhoo(i,j+1)))*(G(i,j+1)-G(i,j)))-...
                    TW_o(i,j)*(Po(i,j)-Po(i,j-1)-1/144*(1/2*(rhoo(i,j)+rhoo(i,j-1)))*(G(i,j)-G(i,j-1)))+Qo(i,j)-...
                    dx(i,j)*dy(i,j)*dz(i,j)/5.615/dt*phi_new(i,j)*(1-Sw_new(i,j)-Sg_new(i,j))/Bo(i,j)+dx(i,j)*dy(i,j)*dz(i,j)/5.615/dt*phi_old(i,j)*(1-Sw_old(i,j)-Sg_old(i,j))/Bo_old(i,j);
            
            Rw(i,j)=TN_w(i,j)*(Pw(i+1,j)-Pw(i,j)-1/144*(0.5*(rhow(i,j)+rhow(i+1,j)))*(G(i+1,j)-G(i,j)))-...
                    TS_w(i,j)*(Pw(i,j)-Pw(i-1,j)-1/144*(0.5*(rhow(i,j)+rhow(i-1,j)))*(G(i,j)-G(i-1,j)))+...
                    TE_w(i,j)*(Pw(i,j+1)-Pw(i,j)-1/144*(0.5*(rhow(i,j)+rhow(i,j+1)))*(G(i,j+1)-G(i,j)))-...
                    TW_w(i,j)*(Pw(i,j)-Pw(i,j-1)-1/144*(0.5*(rhow(i,j)+rhow(i,j-1)))*(G(i,j)-G(i,j-1)))+Qw(i,j)-...
                    (dx(i,j)*dy(i,j)*dz(i,j))/5.615/dt*phi_new(i,j)*Sw_new(i,j)/Bw(i,j)+(dx(i,j)*dy(i,j)*dz(i,j))/5.615/dt*phi_old(i,j)*Sw_old(i,j)/Bw_old(i,j);
            
            Rg(i,j)=TN_g(i,j)*(Pg(i+1,j)-Pg(i,j)-1/144*(1/2*(rhog(i,j)+rhog(i+1,j)))*(G(i+1,j)-G(i,j)))-...
                    TS_g(i,j)*(Pg(i,j)-Pg(i-1,j)-1/144*(1/2*(rhog(i,j)+rhog(i-1,j)))*(G(i,j)-G(i-1,j)))+...
                    TE_g(i,j)*(Pg(i,j+1)-Pg(i,j)-1/144*(1/2*(rhog(i,j)+rhog(i,j+1)))*(G(i,j+1)-G(i,j)))-...
                    TW_g(i,j)*(Pg(i,j)-Pg(i,j-1)-1/144*(1/2*(rhog(i,j)+rhog(i,j-1)))*(G(i,j)-G(i,j-1)))+...
                    TN_o(i,j)*(1/2*(Rso(i+1,j)+Rso(i,j)))*(Po(i+1,j)-Po(i,j)-1/144*(1/2*(rhoo(i,j)+rhoo(i+1,j)))*(G(i+1,j)-G(i,j)))-...
                    TS_o(i,j)*(1/2*(Rso(i,j)+Rso(i-1,j)))*(Po(i,j)-Po(i-1,j)-1/144*(1/2*(rhoo(i,j)+rhoo(i-1,j)))*(G(i,j)-G(i-1,j)))+...
                    TE_o(i,j)*(1/2*(Rso(i,j+1)+Rso(i,j)))*(Po(i,j+1)-Po(i,j)-1/144*(1/2*(rhoo(i,j)+rhoo(i,j+1)))*(G(i,j+1)-G(i,j)))-...
                    TW_o(i,j)*(1/2*(Rso(i,j)+Rso(i,j-1)))*(Po(i,j)-Po(i,j-1)-1/144*(1/2*(rhoo(i,j)+rhoo(i,j-1)))*(G(i,j)-G(i,j-1)))+Qg(i,j)-...
                    (dx(i,j)*dy(i,j)*dz(i,j)/5.615/dt*(phi_new(i,j)*Sg_new(i,j)/Bg(i,j)-phi_old(i,j)*Sg_old(i,j)/Bg_old(i,j)+...
                                                       phi_new(i,j)*Rso(i,j)/Bo(i,j)*(1-Sw_new(i,j)-Sg_new(i,j))-phi_old(i,j)*Rso_old(i,j)/Bo_old(i,j)*(1-Sw_old(i,j)-Sg_old(i,j))));     
           
        end
    end
end

end