function [rhov,Bov,muov,Rsov]=PVT_oil(Pv,PVT_OIL)
% This function is used linear extroplation for the PVT data for oil based
% on the PVT data listed in Table 1
% !This subroutine only apply to the data that follows the normal sequence
P=PVT_OIL(:,1);
rho=PVT_OIL(:,2);
Bo=PVT_OIL(:,3);
muo=PVT_OIL(:,4);
Rso=PVT_OIL(:,5);

% Linear extrapolation
num=1;

while num<=length(P)-1
    if Pv>=P(num) && Pv<P(num+1)
        rhov=rho(num)+(rho(num+1)-rho(num))*(Pv-P(num))/(P(num+1)-P(num));
        Bov=Bo(num)+(Bo(num+1)-Bo(num))*(Pv-P(num))/(P(num+1)-P(num));
        muov=muo(num)+(muo(num+1)-muo(num))*(Pv-P(num))/(P(num+1)-P(num));
        Rsov=Rso(num)+(Rso(num+1)-Rso(num))*(Pv-P(num))/(P(num+1)-P(num));
    end
    num=num+1;
end

end
