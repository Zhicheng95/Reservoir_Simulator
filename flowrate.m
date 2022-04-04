function [Qo,Qw,Qg,Qt,Psf]=flowrate(sp_condition,sp_data,PVT_WATER,PVT_OIL,PVT_GAS,OW,OG,numi,numj,Po,Sw,Sg,index,kx,ky,dx,dy,dz)
% This function is used calcualte the flow rates of the wells based on the
% different conditions. 
% Different conditions are listed in the switch function. In total 8
% conditions are listed. 

% In the input variables, sp_condition is the matrix of (numi,numj) which
% contains the elements rage from 1 to 8 and represents the following
% senarios. 

% sp_data contains the elements which reflects the data that should be
% given to that grid block. And sp_condition and sp_data these two matrixes
% correpond to each other.


Bw=zeros(numi,numj);
muw=zeros(numi,numj);
Bo=zeros(numi,numj);
muo=zeros(numi,numj);
Rso=zeros(numi,numj);
Bg=zeros(numi,numj);
mug=zeros(numi,numj);
krw=zeros(numi,numj);
krg=zeros(numi,numj);
kro=zeros(numi,numj);
Pcow=zeros(numi,numj);
Pcgo=zeros(numi,numj);

Qo=zeros(numi,numj);
Qg=zeros(numi,numj);
Qw=zeros(numi,numj);
Qt=zeros(numi,numj);
Qliq=zeros(numi,numj);
Psf=zeros(numi,numj);
% Oil flow rate is in STB/D
% Water flow/injection rate is in STB/D
% Gas flow/injection rate is in SCF/D
% Total flow rate is in STB/D
% Bottomhole sandface pressure is in PSIA
alphac=5.615;
betac=-1.127*10^(-3);

 %1. When total flow rate is specified
 %2. When oil flow rate is specified
 %3. When water flow rate is specified
 %4. When gas flow rate is specified
 %5. When total liquid flow rate is specified
 %6. When Water injection flow rate is specified
 %7. When Gas injection flow rate is specified
 %8. When Psf is specified
 
for i=1:numi
    for j=1:numj
        if index(i,j)==2
           [krw(i,j),krg(i,j),kro(i,j),Pcow(i,j),Pcgo(i,j)]=relaperm(Sw(i,j),Sg(i,j),OW,OG);
           [~,Bw(i,j),muw(i,j)]=PVT_water(Po(i,j)-Pcow(i,j),PVT_WATER);
           [~,Bo(i,j),muo(i,j),Rso(i,j)]=PVT_oil(Po(i,j),PVT_OIL);
           [~,Bg(i,j),mug(i,j)]=PVT_gas(Po(i,j)+Pcgo(i,j),PVT_GAS);
        end
    end
end

for i=1:numi
    for j=1:numj
        if index(i,j)==2
           
        switch sp_condition(i,j)
            case 1 % When Q total is specified
                Qt(i,j)=sp_data(i,j);
                Qo(i,j)=(Qt(i,j)*kro(i,j)/muo(i,j)/Bo(i,j))/(kro(i,j)/muo(i,j)/Bo(i,j)+...
                    krw(i,j)/muw(i,j)/Bw(i,j)+krg(i,j)/mug(i,j)/Bg(i,j)/alphac+Rso(i,j)*kro(i,j)/alphac/muo(i,j)/Bo(i,j));
                
                Qw(i,j)=(Qt(i,j)*krw(i,j)/muw(i,j)/Bw(i,j))/(kro(i,j)/muo(i,j)/Bo(i,j)+...
                    krw(i,j)/muw(i,j)/Bw(i,j)+krg(i,j)/mug(i,j)/Bg(i,j)/alphac+Rso(i,j)*kro(i,j)/alphac/muo(i,j)/Bo(i,j));
                
                Qg(i,j)=Qt(i,j)*(krg(i,j)/mug(i,j)/Bg(i,j)+Rso(i,j)*alphac*kro(i,j)/muo(i,j)/Bo(i,j))/(kro(i,j)/muo(i,j)/Bo(i,j)+...
                    krw(i,j)/muw(i,j)/Bw(i,j)+krg(i,j)/mug(i,j)/Bg(i,j)/alphac+Rso(i,j)*kro(i,j)/alphac/muo(i,j)/Bo(i,j));
                
            case 2 % When Qo is specified 
                Qo(i,j)=sp_data(i,j);
                
                Qw(i,j)=Qo(i,j)*krw(i,j)/muw(i,j)/Bw(i,j)/(kro(i,j)/muo(i,j)/Bo(i,j));
                
                Qg(i,j)=Qo(i,j)*alphac*(krg(i,j)/mug(i,j)/Bg(i,j)/alphac+kro(i,j)*Rso(i,j)/muo(i,j)/Bo(i,j)/alphac)/(kro(i,j)/muo(i,j)/Bo(i,j));

            case 3 % When Qw si specified
                Qw(i,j)=sp_data(i,j);
                
                Qo(i,j)=Qw(i,j)*(kro(i,j)/muo(i,j)/Bo(i,j)/(krw(i,j)/muw(i,j)/Bw(i,j)));
                
                Qg(i,j)=Qw(i,j)*(krg(i,j)/mug(i,j)/Bg(i,j)+Rso(i,j)*kro(i,j)/mug(i,j)/Bg(i,j))/(krw(i,j)/muw(i,j)/Bw(i,j));
            case 4 % When Qg is specified
                Qg(i,j)=sp_data(i,j);
                
                Qo(i,j)=Qg(i,j)*(kro(i,j)/muo(i,j)/Bo(i,j))/(krg(i,j)/mug(i,j)/Bg(i,j)+Rso(i,j)*kro(i,j)/muo(i,j)/Bo(i,j));
                
                Qw(i,j)=Qg(i,j)*(krw(i,j)/muw(i,j)/Bw(i,j))/(krg(i,j)/mug(i,j)/Bg(i,j)+Rso(i,j)*kro(i,j)/muo(i,j)/Bo(i,j));
            case 5 % When total Qliq is specified
                Qliq(i,j)=sp_data(i,j);
                
                Qo(i,j)=kro(i,j)/muo(i,j)/Bo(i,j)/(kro(i,j)/muo(i,j)/Bo(i,j)+krw(i,j)/muw(i,j)/Bw(i,j))*Qliq(i,j);
                
                Qw(i,j)=krw(i,j)/muw(i,j)/Bw(i,j)/(kro(i,j)/muo(i,j)/Bo(i,j)+krw(i,j)/muw(i,j)/Bw(i,j))*Qliq(i,j);
                
                Qg(i,j)=(krg(i,j)/mug(i,j)/Bg(i,j)+Rso(i,j)*kro(i,j)/muo(i,j)/Bo(i,j))/(kro(i,j)/muo(i,j)/Bo(i,j)+krw(i,j)/muw(i,j)/Bw(i,j))*Qliq(i,j);
            case 6 % When Water injection flow rate is specified
                Qw(i,j)=sp_data(i,j);
                
                Qo(i,j)=0;
                
                Qg(i,j)=0;
            case 7 % When Gas injection flow rate is specified
                Qg(i,j)=sp_data(i,j);
                
                Qo(i,j)=0;
                
                Qw(i,j)=0;
            case 8 % When Psf is specified
                Psf(i,j)=sp_data(i,j);
                
                rw=0.25; % well radius is 0.25 ft     
                
                re=0.28*(((ky(i,j)/kx(i,j))^(1/2))*(dx(i,j)^2)+((kx(i,j)/ky(i,j))^(1/2))*(dy(i,j)^2))/((ky(i,j)/kx(i,j)^(1/4)+(kx(i,j)/ky(i,j))^(1/4))); 
               
                
                Qo(i,j)=betac*((2*pi*((kx(i,j)*ky(i,j))^(1/2))*dz(i,j)*kro(i,j))/(muo(i,j)*Bo(i,j)*log(re/rw)))*(Po(i,j)-Psf(i,j));
                
                Qw(i,j)=betac*((2*pi*((kx(i,j)*ky(i,j))^(1/2))*dz(i,j)*krw(i,j))/(muw(i,j)*Bw(i,j)*log(re/rw)))*(Po(i,j)-Pcow(i,j)-Psf(i,j)); 
                
                Qg(i,j)=betac*((2*pi*((kx(i,j)*ky(i,j))^(1/2))*dz(i,j)*krw(i,j))/(mug(i,j)*Bg(i,j)*log(re/rw)))*(Po(i,j)+Pcgo(i,j)-Psf(i,j))+ ...
                            betac*Rso(i,j)*((2*pi*((kx(i,j)*ky(i,j))^(1/2))*dz(i,j)*kro(i,j))/(muo(i,j)*Bo(i,j)*log(re/rw)))*(Po(i,j)-Psf(i,j));    
                
        end
        end
    end
end
end
