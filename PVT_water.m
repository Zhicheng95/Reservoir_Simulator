function [rhowv,Bwv,muwv]=PVT_water(Pv,PVT_WATER)
% This function is used linear extroplation for the PVT data for oil based
% on the PVT data listed in Table 1
% !This subroutine only apply to the data that follows the normal sequence
P=PVT_WATER(:,1);
rhow=PVT_WATER(:,2);
Bw=PVT_WATER(:,3);
muw=PVT_WATER(:,4);

% Linear extrapolation
num=1;

while num<=length(P)-1
    if Pv>=P(num) && Pv<P(num+1)
        rhowv=rhow(num)+(rhow(num+1)-rhow(num))*(Pv-P(num))/(P(num+1)-P(num));
        Bwv=Bw(num)+(Bw(num+1)-Bw(num))*(Pv-P(num))/(P(num+1)-P(num));
        muwv=muw(num)+(muw(num+1)-muw(num))*(Pv-P(num))/(P(num+1)-P(num));
       
    end
    num=num+1;
end

end