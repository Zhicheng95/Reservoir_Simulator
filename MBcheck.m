function [tag_inc,tag_com,oil_inc,water_inc,gas_inc,oil_com,water_com,gas_com]=MBcheck(numi,numj,numt,time,Po,Po_old,Sw_old,Sg_old,Po_ini,phi_ini,...
    Sw_ini,Sg_ini,Sw_new,Sg_new,Qo,Qw,Qg,Qov,Qwv,Qgv,PVT_GAS,PVT_WATER,PVT_OIL,OW,OG,CR,dt,index,dx,dy,dz,Bo_ini,Bw_ini,Bg_ini)

Bo=zeros(numi,numj);
Rso=zeros(numi,numj);
Pcow=zeros(numi,numj);
Pcow_old=zeros(numi,numj);
Pcgo=zeros(numi,numj);
Pcgo_old=zeros(numi,numj);
Bw=zeros(numi,numj);
Bg=zeros(numi,numj);
Bo_old=zeros(numi,numj);
Bw_old=zeros(numi,numj);
Bg_old=zeros(numi,numj);
Rso_old=zeros(numi,numj);
phi_new=zeros(numi,numj);
phi_old=zeros(numi,numj);

Vb=zeros(numi,numj);
So_new=zeros(numi,numj);
So_old=zeros(numi,numj);
So_ini=zeros(numi,numj);


MB_inc_oil_new=zeros(numi,numj);
MB_inc_oil_old=zeros(numi,numj);
MB_inc_gas_new=zeros(numi,numj);
MB_inc_gas_old=zeros(numi,numj);
MB_inc_water_new=zeros(numi,numj);
MB_inc_water_old=zeros(numi,numj);
MB_well_oil=zeros(numi,numj);
MB_well_gas=zeros(numi,numj);
MB_well_water=zeros(numi,numj);

for i=1:numi
    for j=1:numj
        if index(i,j)~=0
            [~,Bo(i,j),~,Rso(i,j)]=PVT_oil(Po(i,j),PVT_OIL);
            [~,~,~,Pcow(i,j),Pcgo(i,j)]=relaperm(Sw_new(i,j),Sg_new(i,j),OW,OG);
            [~,Bw(i,j),~]=PVT_water(Po(i,j)-Pcow(i,j),PVT_WATER);
            [~,Bg(i,j),~]=PVT_gas(Po(i,j)+Pcgo(i,j),PVT_GAS);
            
            
            [~,~,~,Pcow_old(i,j),Pcgo_old(i,j)]=relaperm(Sw_old(i,j),Sg_old(i,j),OW,OG);
            [~,Bo_old(i,j),~,Rso_old(i,j)]=PVT_oil(Po_old(i,j),PVT_OIL);
            [~,Bw_old(i,j),~]=PVT_water(Po_old(i,j)-Pcow_old(i,j),PVT_WATER);
            [~,Bg_old(i,j),~]=PVT_gas(Po_old(i,j)+Pcgo_old(i,j),PVT_GAS);
            phi_new(i,j)=phi_ini(i,j)*(1+CR*(Po(i,j)-Po_ini(i,j)));
            phi_old(i,j)=phi_ini(i,j)*(1+CR*(Po_old(i,j)-Po_ini(i,j)));
        end
    end
end


%% Incremental MBC

for i=1:numi
    for j=1:numj
        if index(i,j)~=0
            Vb(i,j)=dx(i,j)*dy(i,j)*dz(i,j);
            So_new(i,j)=1-Sw_new(i,j)-Sg_new(i,j);
            So_old(i,j)=1-Sw_old(i,j)-Sg_old(i,j);
            
            MB_inc_oil_new(i,j)=Vb(i,j)*phi_new(i,j)*So_new(i,j)/Bo(i,j)/5.615;
            MB_inc_oil_old(i,j)=Vb(i,j)*phi_old(i,j)*So_old(i,j)/Bo_old(i,j)/5.615;
            
            MB_inc_gas_new(i,j)=Vb(i,j)*phi_new(i,j)*Sg_new(i,j)/Bg(i,j)/5.615;
            MB_inc_gas_old(i,j)=Vb(i,j)*phi_old(i,j)*Sg_old(i,j)/Bg_old(i,j)/5.615;
            
            MB_inc_water_new(i,j)=Vb(i,j)*phi_new(i,j)*Sw_new(i,j)/Bw(i,j)/5.615;
            MB_inc_water_old(i,j)=Vb(i,j)*phi_old(i,j)*Sw_old(i,j)/Bw_old(i,j)/5.615;
                if index(i,j)==2
                    MB_well_oil(i,j)=Qo(i,j)*dt;
                    MB_well_gas(i,j)=Qg(i,j)*dt;
                    MB_well_water(i,j)=Qw(i,j)*dt;
                end
        end
    end
end

tag_oil_inc=0;
tag_gas_inc=0;
tag_water_inc=0;
tag_inc=0;

oil_inc=(sum(sum(MB_inc_oil_new))-sum(sum(MB_inc_oil_old)))/sum(sum(MB_well_oil));
if abs(oil_inc-1)<10^(-5)
    tag_oil_inc=1;
end
water_inc=(sum(sum(MB_inc_water_new))-sum(sum(MB_inc_water_old)))/sum(sum(MB_well_water));
if abs(water_inc-1)<10^(-5)
    tag_water_inc=1;
end

gas_inc=(sum(sum(MB_inc_oil_new))-sum(sum(MB_inc_oil_old)))/sum(sum(MB_well_oil));
if abs(gas_inc-1)<10^(-5)
    tag_gas_inc=1;
end

if tag_oil_inc==1 && tag_water_inc==1 && tag_gas_inc==1
    tag_inc=1;
end

%% Comulative MBC

MB_com_oil_ini=zeros(numi,numj);
MB_com_gas_ini=zeros(numi,numj);
MB_com_water_ini=zeros(numi,numj);

MB_com_oil=zeros(numi,numj,numt+1);
MB_com_gas=zeros(numi,numj,numt+1);
MB_com_water=zeros(numi,numj,numt+1);

MB_com_well_oil=zeros(numi,numj,numt+1);
MB_com_well_gas=zeros(numi,numj,numt+1);
MB_com_well_water=zeros(numi,numj,numt+1);

for i=1:numi
    for j=1:numj
        if index(i,j)~=0
            So_ini(i,j)=1-Sg_ini(i,j)-Sw_ini(i,j);
            MB_com_oil_ini(i,j)=Vb(i,j)*phi_ini(i,j)*So_ini(i,j)/Bo_ini(i,j)/5.615;                        
            MB_com_gas_ini(i,j)=Vb(i,j)*phi_ini(i,j)*Sg_ini(i,j)/Bg_ini(i,j)/5.615;
            MB_com_water_ini(i,j)=Vb(i,j)*phi_ini(i,j)*Sw_ini(i,j)/Bw_ini(i,j)/5.615;
            
%             MB_com_oil(i,j,1)=MB_com_oil_ini(i,j);       
%             MB_com_gas(i,j,1)=MB_com_gas_ini(i,j);
%             MB_com_water(i,j,1)=MB_com_water_ini(i,j);
%                     if index(i,j)==2
%                         MB_com_well_oil(i,j,1)=Qov(i,j,1)*dt;
%                         MB_com_well_gas(i,j,1)=Qgv(i,j,1)*dt;
%                         MB_com_well_water(i,j,1)=Qwv(i,j,1)*dt;
%                     end
         end
    end
end
            MB_com_oil(:,:,1)=MB_com_oil_ini;       
            MB_com_gas(:,:,1)=MB_com_gas_ini;
            MB_com_water(:,:,1)=MB_com_water_ini;

for t=2:time
    for i=1:numi
      for j=1:numj
              if index(i,j)~=0
                  
                Vb(i,j)=dx(i,j)*dy(i,j)*dz(i,j);
                So_new(i,j)=1-Sw_new(i,j)-Sg_new(i,j);          
                                       
                MB_com_oil(i,j,t)=Vb(i,j)*phi_new(i,j)*So_new(i,j)/Bo(i,j)/5.615;
                       
                MB_com_gas(i,j,t)=Vb(i,j)*phi_new(i,j)*Sg_new(i,j)/Bg(i,j)/5.615;
                
                MB_com_water(i,j,t)=Vb(i,j)*phi_new(i,j)*Sw_new(i,j)/Bw(i,j)/5.615;
                    if index(i,j)==2
                       
                        MB_com_well_oil(i,j,t)=Qov(i,j,t)*dt;
                        MB_com_well_gas(i,j,t)=Qgv(i,j,t)*dt;
                        MB_com_well_water(i,j,t)=Qwv(i,j,t)*dt;
                    end
              end
       end
     end
end

tag_oil_com=0;
tag_gas_com=0;
tag_water_com=0;
tag_com=0;

oil_com=(sum(sum(MB_com_oil(:,:,time)))-sum(sum(MB_com_oil_ini)))/sum(sum(sum(MB_com_well_oil)));
if abs(oil_com-1)<10^(-5)
    tag_oil_com=1;
end

water_com=(sum(sum(MB_com_water(:,:,time)))-sum(sum(MB_com_water_ini)))/sum(sum(sum(MB_com_well_water)));
if abs(water_com-1)<10^(-5)
    tag_water_com=1;
end

gas_com=(sum(sum(MB_com_oil(:,:,time)))-sum(sum(MB_com_oil_ini)))/sum(sum(sum(MB_com_well_oil)));
if abs(gas_com-1)<10^(-5)
    tag_gas_com=1;
end

if (tag_oil_com==1 && tag_water_com==1 && tag_gas_com==1)
    tag_com=1;
end
end





        
        