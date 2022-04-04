function [rhogv,Bgv,mugv]=PVT_gas(Pv,PVT_GAS)
% This function is used linear extroplation for the PVT data for oil based
% on the PVT data listed in Table 1
% !This subroutine only apply to the data that follows the normal sequence
P=PVT_GAS(:,1);
rhog=PVT_GAS(:,2);
Bg=PVT_GAS(:,3);
mug=PVT_GAS(:,4);

% Linear extrapolation
num=1;

while num<=length(P)-1
    if Pv>=P(num) && Pv<P(num+1)
        rhogv=rhog(num)+(rhog(num+1)-rhog(num))*(Pv-P(num))/(P(num+1)-P(num));
        Bgv=Bg(num)+(Bg(num+1)-Bg(num))*(Pv-P(num))/(P(num+1)-P(num));
        mugv=mug(num)+(mug(num+1)-mug(num))*(Pv-P(num))/(P(num+1)-P(num));
       
    end
    num=num+1;
end
end