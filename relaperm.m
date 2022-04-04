function [krwv,krgv,krov,Pcowv,Pcgov]=relaperm(Swv,Sgv,OW,OG)
% This function is used to calculate the relative permeability of 3 phases
% system for each grid block based on the Stone's second model.
% The input is the water % gas saturation and the output is the relative
% permeabilities and the capilary pressures.
Sw=OW(:,1);
krw=OW(:,2);
krow=OW(:,3);
Pcow=OW(:,4);
Swirr=min(Sw);

Sg=OG(:,1);
krg=OG(:,2);
krog=OG(:,3);
Pcgo=OG(:,4);
Sor=1-max(Sg);

krwv=0;
krgv=0;
krowv=0;
krogv=0;

if Swv<=Swirr
    krwv=0;
    krowv=1;
end

if Sgv>=1-Sor
    krgv=1;
    krogv=0;
end

numw=1;
numg=1;
while numw<=length(Sw)-1
    if Swv>=Sw(numw) && Swv<Sw(numw+1)
        krwv=krw(numw)+(krw(numw+1)-krw(numw))*(Swv-Sw(numw))/(Sw(numw+1)-Sw(numw));
        krowv=krow(numw)+(krow(numw+1)-krow(numw))*(Swv-Sw(numw))/(Sw(numw+1)-Sw(numw));
        Pcowv=Pcow(numw)+(Pcow(numw+1)-Pcow(numw))*(Swv-Sw(numw))/(Sw(numw+1)-Sw(numw));
    end
    numw=numw+1;
end

while numg<=length(Sg)-1
    if Sgv>=Sg(numg) && Sgv<Sg(numg+1)
        krgv=krg(numg)+(krg(numg+1)-krg(numg))*(Sgv-Sg(numg))/(Sg(numg+1)-Sg(numg));
        krogv=krog(numg)+(krog(numg+1)-krog(numg))*(Sgv-Sg(numg))/(Sg(numg+1)-Sg(numg));
        Pcgov=Pcgo(numg)+(Pcgo(numg+1)-Pcgo(numg))*(Sgv-Sg(numg))/(Sg(numg+1)-Sg(numg));
    end
    numg=numg+1;
end

krov=(krowv+krwv)*(krogv+krgv)-(krwv+krgv);
if krov<0
    krov=0;
end
end

