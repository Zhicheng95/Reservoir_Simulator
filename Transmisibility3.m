function [TN_o,TN_w,TN_g,TS_o,TS_w,TS_g,TE_o,TE_w,TE_g,TW_o,TW_w,TW_g]=Transmisibility3(Po,numi,numj,index,G,dx,dy,dz,kx,ky,rhoo,Bo,muo,...
   krw,krg,kro,Pcow,Pcgo,rhow,Bw,muw,rhog,Bg,mug)
% This function is used to calculate the transmisibility of each grid block
% The unit of permeability is mD.
betac=1.127e-3;
% The input for this function is Sw, Sg and Pressure for the whole
% reservoir
Pot_o=zeros(numi,numj);
Pot_g=zeros(numi,numj);
Pot_w=zeros(numi,numj);

TN_o=zeros(numi,numj);
TN_w=zeros(numi,numj);
TN_g=zeros(numi,numj);

TS_o=zeros(numi,numj);
TS_w=zeros(numi,numj);
TS_g=zeros(numi,numj);

TE_o=zeros(numi,numj);
TE_w=zeros(numi,numj);
TE_g=zeros(numi,numj);

TW_o=zeros(numi,numj);
TW_w=zeros(numi,numj);
TW_g=zeros(numi,numj);

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

AKDN=zeros(numi,numj);
AKDS=zeros(numi,numj);
AKDE=zeros(numi,numj);
AKDW=zeros(numi,numj);

%% For All the transmisibilities inside    
for i=1:numi % numx,numy are the number of blocks which are considered as the real reservoir
    for j=1:numj
        if i==1 || i==numi || j==1 || j==numj
            Pot_o(i,j)=0;
            Pot_g(i,j)=0;
            Pot_w(i,j)=0;
        end
    end
end

for i=1:numi
    for j=1:numj
        
        if index(i,j)~=0
%             [rhoo(i,j),Bo(i,j),muo(i,j),Rso(i,j)]=PVT_oil(Po(i,j),PVT_OIL);
%             [krw(i,j),krg(i,j),kro(i,j),Pcow(i,j),Pcgo(i,j)]=relaperm(Sw(i,j),Sg(i,j),OW,OG);
%             [rhow(i,j),Bw(i,j),muw(i,j)]=PVT_water(Po(i,j)-Pcow(i,j),PVT_WATER);
%             [rhog(i,j),Bg(i,j),mug(i,j)]=PVT_gas(Po(i,j)+Pcgo(i,j),PVT_GAS);
            
            Pot_o(i,j)=Po(i,j)-1/144*rhoo(i,j)*G(i,j);
            Pot_g(i,j)=Po(i,j)+Pcgo(i,j)-1/144*rhog(i,j)*G(i,j);
            Pot_w(i,j)=Po(i,j)-Pcow(i,j)-1/144*rhow(i,j)*G(i,j);
        
            AKDN(i,j)=2/(dx(i,j)/(dy(i,j)*dz(i,j)*kx(i,j))+dx(i+1,j)/(dy(i+1,j)*dz(i+1,j)*kx(i+1,j))); %k_x*A_x/dx, harmonic average
            AKDS(i,j)=2/(dx(i,j)/(dy(i,j)*dz(i,j)*kx(i,j))+dx(i-1,j)/(dy(i-1,j)*dz(i-1,j)*kx(i-1,j)));
            AKDE(i,j)=2/(dy(i,j)/(dx(i,j)*dz(i,j)*ky(i,j))+dy(i,j+1)/(dx(i,j+1)*dz(i,j+1)*ky(i,j+1)));
            AKDW(i,j)=2/(dy(i,j)/(dx(i,j)*dz(i,j)*ky(i,j))+dy(i,j-1)/(dx(i,j-1)*dz(i,j-1)*ky(i,j-1)));
        end 
    end
end

for i=2:numi-1
    for j=2:numj-1
            %%% Oil Transmisibility
            
            if Pot_o(i+1,j)>=Pot_o(i,j)
                TN_o(i,j)=betac*AKDN(i,j)*kro(i+1,j)/(1/2*(muo(i+1,j)+muo(i,j)))/(1/2*(Bo(i+1,j)+Bo(i,j)));
            else
                TN_o(i,j)=betac*AKDN(i,j)*kro(i,j)/(1/2*(muo(i+1,j)+muo(i,j)))/(1/2*(Bo(i+1,j)+Bo(i,j)));
            end
            
            if Pot_o(i-1,j)>=Pot_o(i,j)
                TS_o(i,j)=betac*AKDS(i,j)*kro(i-1,j)/(1/2*(muo(i-1,j)+muo(i,j)))/(1/2*(Bo(i-1,j)+Bo(i,j)));
            else
                TS_o(i,j)=betac*AKDS(i,j)*kro(i,j)/(1/2*(muo(i-1,j)+muo(i,j)))/(1/2*(Bo(i-1,j)+Bo(i,j)));
            end
             
            if Pot_o(i,j+1)>=Pot_o(i,j)
                TE_o(i,j)=betac*AKDE(i,j)*kro(i,j+1)/(1/2*(muo(i,j+1)+muo(i,j)))/(1/2*(Bo(i,j+1)+Bo(i,j)));
            else
                TE_o(i,j)=betac*AKDE(i,j)*kro(i,j)/(1/2*(muo(i,j+1)+muo(i,j)))/(1/2*(Bo(i,j+1)+Bo(i,j)));
            end
            
            if Pot_o(i,j-1)>=Pot_o(i,j)
                TW_o(i,j)=betac*AKDW(i,j)*kro(i,j-1)/(1/2*(muo(i,j-1)+muo(i,j)))/(1/2*(Bo(i,j-1)+Bo(i,j)));
            else
                TW_o(i,j)=betac*AKDW(i,j)*kro(i,j)/(1/2*(muo(i,j-1)+muo(i,j)))/(1/2*(Bo(i,j-1)+Bo(i,j)));
            end
            
            %%% Water Transmisibility
            
            if Pot_w(i+1,j)>=Pot_w(i,j)
                TN_w(i,j)=betac*AKDN(i,j)*krw(i+1,j)/(1/2*(muw(i+1,j)+muw(i,j)))/(1/2*(Bw(i+1,j)+Bw(i,j)));
            else
                TN_w(i,j)=betac*AKDN(i,j)*krw(i,j)/(1/2*(muw(i+1,j)+muw(i,j)))/(1/2*(Bw(i+1,j)+Bw(i,j)));
            end
            
            if Pot_w(i-1,j)>=Pot_w(i,j)
                TS_w(i,j)=betac*AKDS(i,j)*krw(i-1,j)/(1/2*(muw(i-1,j)+muw(i,j)))/(1/2*(Bw(i-1,j)+Bw(i,j)));
            else
                TS_w(i,j)=betac*AKDS(i,j)*krw(i,j)/(1/2*(muw(i-1,j)+muw(i,j)))/(1/2*(Bw(i-1,j)+Bw(i,j)));
            end
             
            if Pot_w(i,j+1)>=Pot_w(i,j)
                TE_w(i,j)=betac*AKDE(i,j)*krw(i,j+1)/(1/2*(muw(i,j+1)+muw(i,j)))/(1/2*(Bw(i,j+1)+Bw(i,j)));
            else
                TE_w(i,j)=betac*AKDE(i,j)*krw(i,j)/(1/2*(muw(i,j+1)+muw(i,j)))/(1/2*(Bw(i,j+1)+Bw(i,j)));
            end
            
            if Pot_w(i,j-1)>=Pot_w(i,j)
                TW_w(i,j)=betac*AKDW(i,j)*krw(i,j-1)/(1/2*(muw(i,j-1)+muw(i,j)))/(1/2*(Bw(i,j-1)+Bw(i,j)));
            else
                TW_w(i,j)=betac*AKDW(i,j)*krw(i,j)/(1/2*(muw(i,j-1)+muw(i,j)))/(1/2*(Bw(i,j-1)+Bw(i,j)));
            end
            
            %%% Gas Transmisibility
            
            if Pot_g(i+1,j)>=Pot_g(i,j)
                TN_g(i,j)=betac*AKDN(i,j)*krg(i+1,j)/(1/2*(mug(i+1,j)+mug(i,j)))/(1/2*(Bg(i+1,j)+Bg(i,j)));
            else
                TN_g(i,j)=betac*AKDN(i,j)*krg(i,j)/(1/2*(mug(i+1,j)+mug(i,j)))/(1/2*(Bg(i+1,j)+Bg(i,j)));
            end
            
            if Pot_g(i-1,j)>=Pot_g(i,j)
                TS_g(i,j)=betac*AKDS(i,j)*krg(i-1,j)/(1/2*(mug(i-1,j)+mug(i,j)))/(1/2*(Bg(i-1,j)+Bg(i,j)));
            else
                TS_g(i,j)=betac*AKDS(i,j)*krg(i,j)/(1/2*(mug(i-1,j)+mug(i,j)))/(1/2*(Bg(i-1,j)+Bg(i,j)));
            end
             
            if Pot_g(i,j+1)>=Pot_g(i,j)
                TE_g(i,j)=betac*AKDE(i,j)*krg(i,j+1)/(1/2*(mug(i,j+1)+mug(i,j)))/(1/2*(Bg(i,j+1)+Bg(i,j)));
            else
                TE_g(i,j)=betac*AKDE(i,j)*krg(i,j)/(1/2*(mug(i,j+1)+mug(i,j)))/(1/2*(Bg(i,j+1)+Bg(i,j)));
            end
            
            if Pot_g(i,j-1)>=Pot_g(i,j)
                TW_g(i,j)=betac*AKDW(i,j)*krg(i,j-1)/(1/2*(mug(i,j-1)+mug(i,j)))/(1/2*(Bg(i,j-1)+Bg(i,j)));
            else
                TW_g(i,j)=betac*AKDW(i,j)*krg(i,j)/(1/2*(mug(i,j-1)+mug(i,j)))/(1/2*(Bg(i,j-1)+Bg(i,j)));
            end
            
    end 

end

%% Boundary 
for j=2:numj-1
    for i=2:numi-1
        if index(i,j)==0
            TS_o(i,j)=0;TN_o(i,j)=0;TE_o(i,j)=0;TW_o(i,j)=0;
            TS_g(i,j)=0;TN_g(i,j)=0;TE_g(i,j)=0;TW_g(i,j)=0; 
            TS_w(i,j)=0;TN_w(i,j)=0;TE_w(i,j)=0;TW_w(i,j)=0;
        end
            
        if index(i,j)~=0 && index(i+1,j)==0 % Which means the north side of i+1,j is boundary
            TN_o(i,j)=0;TN_g(i,j)=0;TN_w(i,j)=0;
        end
        
        if index(i,j)~=0 && index(i-1,j)==0 % Which means the north side of i+1,j is boundary
            TS_o(i,j)=0;TS_g(i,j)=0;TS_w(i,j)=0;
        end
        
        if index(i,j)~=0 && index(i,j+1)==0 % Which means the north side of i+1,j is boundary
            TE_o(i,j)=0;TE_g(i,j)=0;TE_w(i,j)=0;
        end
        
        if index(i,j)~=0 && index(i,j-1)==0 % Which means the north side of i+1,j is boundary
            TW_o(i,j)=0;TW_g(i,j)=0;TW_w(i,j)=0;
        end
    end
end

        
% for j=2:numj-1
%     for i=1:numi-1
%         if index(i,j)==0 && index(i+1,j)==1 % Which means the south side of i+1,j is boundary
%             TS_o(i+1,j)=0;
%             TS_w(i+1,j)=0;
%             TS_g(i+1,j)=0;
%         elseif index(i,j)==1 && index(i+1,j)==0 % Which means the north side of i,j is boundary
%             TN_o(i,j)=0;
%             TN_w(i,j)=0;
%             TN_g(i,j)=0;
%         end
%     end
% end
% 
% for i=2:numi-1
%     for j=1:numj-1
%         if index(i,j)==0 && index(i,j+1)==1 % Which means the west side of i,j+1 is boundary
%             TW_o(i,j+1)=0;
%             TW_w(i,j+1)=0;
%             TW_g(i,j+1)=0;
%         elseif index(i,j)==1 && index(i,j+1)==0 % Which means the east side of i,j is boundary
%             TE_o(i,j)=0;
%             TE_w(i,j)=0;
%             TE_g(i,j)=0;
%         end
%     end
% end           

end