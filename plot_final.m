
for i=1:73
    QG1(i)=Qgv(46,22,i);
    QO1(i)=Qov(46,22,i);
    QW1(i)=Qwv(46,22,i);
    PV1(i)=Pov(46,22,i);
    QG2(i)=Qgv(42,16,i);
    QO2(i)=Qov(42,16,i);
    QW2(i)=Qwv(42,16,i);
    PV2(i)=Pov(42,16,i);
    QG3(i)=Qgv(40,10,i);
    QO3(i)=Qov(40,10,i);
    QW3(i)=Qwv(40,10,i);
    PV3(i)=Pov(40,10,i);
    QG4(i)=Qgv(37,23,i);
    QO4(i)=Qov(37,23,i);
    QW4(i)=Qwv(37,23,i);
    PV4(i)=Pov(37,23,i);
    QG5(i)=Qgv(34,12,i);
    QO5(i)=Qov(34,12,i);
    QW5(i)=Qwv(34,12,i);
    PV5(i)=Pov(34,12,i);
    QG6(i)=Qgv(31,15,i);
    QO6(i)=Qov(31,15,i);
    QW6(i)=Qwv(31,15,i);
    PV6(i)=Pov(31,15,i);
    QG7(i)=Qgv(31,12,i);
    QO7(i)=Qov(31,12,i);
    QW7(i)=Qwv(31,12,i);
    PV7(i)=Pov(31,12,i);
    QG8(i)=Qgv(30,6,i);
    QO8(i)=Qov(30,6,i);
    QW8(i)=Qwv(30,6,i);
    PV8(i)=Pov(30,6,i);
    QG9(i)=Qgv(26,19,i);
    QO9(i)=Qov(26,19,i);
    QW9(i)=Qwv(26,19,i);
    PV9(i)=Pov(26,19,i);
    QG10(i)=Qgv(26,7,i);
    QO10(i)=Qov(26,7,i);
    QW10(i)=Qwv(26,7,i);
    PV10(i)=Pov(26,7,i);

end
QWC1(1)=QW1(1);
QOC1(1)=QO1(1);
QGC1(1)=QG1(1);

QWC2(1)=QW2(1);
QOC2(1)=QO2(1);
QGC2(1)=QG2(1);
QWC3(1)=QW3(1);
QOC3(1)=QO3(1);
QGC3(1)=QG3(1);
QWC4(1)=QW4(1);
QOC4(1)=QO4(1);
QGC4(1)=QG4(1);
QWC5(1)=QW5(1);
QOC5(1)=QO5(1);
QGC5(1)=QG5(1);
QWC6(1)=QW6(1);
QOC6(1)=QO6(1);
QGC6(1)=QG6(1);
QWC7(1)=QW7(1);
QOC7(1)=QO7(1);
QGC7(1)=QG7(1);
QWC8(1)=QW8(1);
QOC8(1)=QO8(1);
QGC8(1)=QG8(1);
QWC9(1)=QW9(1);
QOC9(1)=QO9(1);
QGC9(1)=QG9(1);
QWC10(1)=QW10(1);
QOC10(1)=QO10(1);
QGC10(1)=QG10(1);

for i=2:73
         QWC1(i)=QW1(i)+QWC1(i-1);
         QOC1(i)=QO1(i)+QOC1(i-1);
         QGC1(i)=QG1(i)+QGC1(i-1);
         QWC2(i)=QW2(i)+QWC2(i-1);
         QOC2(i)=QO2(i)+QOC2(i-1);
         QGC2(i)=QG2(i)+QGC2(i-1);
         QWC3(i)=QW3(i)+QWC3(i-1);
         QOC3(i)=QO3(i)+QOC3(i-1);
         QGC3(i)=QG3(i)+QGC3(i-1);
         QWC4(i)=QW4(i)+QWC4(i-1);
         QOC4(i)=QO4(i)+QOC4(i-1);
         QGC4(i)=QG4(i)+QGC4(i-1);
         QWC5(i)=QW5(i)+QWC5(i-1);
         QOC5(i)=QO5(i)+QOC5(i-1);
         QGC5(i)=QG5(i)+QGC5(i-1);
         QWC6(i)=QW6(i)+QWC6(i-1);
         QOC6(i)=QO6(i)+QOC6(i-1);
         QGC6(i)=QG6(i)+QGC6(i-1);  
         QWC7(i)=QW7(i)+QWC7(i-1);
         QOC7(i)=QO7(i)+QOC7(i-1);
         QGC7(i)=QG7(i)+QGC7(i-1);
         QWC8(i)=QW8(i)+QWC8(i-1);
         QOC8(i)=QO8(i)+QOC8(i-1);
         QGC8(i)=QG8(i)+QGC8(i-1);
         QWC9(i)=QW9(i)+QWC9(i-1);
         QOC9(i)=QO9(i)+QOC9(i-1);
         QGC9(i)=QG9(i)+QGC9(i-1);
         QWC10(i)=QW5(i)+QWC5(i-1);
         QOC10(i)=QO10(i)+QOC10(i-1);
         QGC10(i)=QG10(i)+QGC10(i-1);
         

end

T=linspace(0,73,73);
T=T*5;
figure(1)
hold on
plot(T,QO1);
plot(T,QO2);
plot(T,QO3);
plot(T,QO4);
plot(T,QO5);
plot(T,QO6);
plot(T,QO7);
plot(T,QO8);
plot(T,QO9);
plot(T,QO10);

hold off
grid on
ylabel('Oil flow rate,STBD')
xlabel('time,day')
title('Fig.1 Oil flow rate vs time')
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')


figure(2)
hold on
plot(T,QG1);
plot(T,QG2);
plot(T,QG3);
plot(T,QG4);
plot(T,QG5);
plot(T,QG6);
plot(T,QG7);
plot(T,QG8);
plot(T,QG9);
plot(T,QG10);
hold off
grid on
ylabel('Gas flow rate,SCFD')
xlabel('time,day')
title('Fig.2 Gas flow rate vs time')
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')

figure(3)
hold on
plot(T,QW1);
plot(T,QW2);
plot(T,QW3);
plot(T,QW4);
plot(T,QW5);
plot(T,QW6);
plot(T,QW7);
plot(T,QW8);
plot(T,QW9);
plot(T,QW10);
hold off
grid on
ylabel('Water flow rate,STBD')
xlabel('time,day')
title('Fig.3 Water flow rate vs time')
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')

figure(4)
hold on
plot(T,QWC1);
plot(T,QWC2);
plot(T,QWC3);
plot(T,QWC4);
plot(T,QWC5);
plot(T,QWC6);
plot(T,QWC7);
plot(T,QWC8);
plot(T,QWC9);
plot(T,QWC10);
hold off
grid on
ylabel('Accumulative Water Production,STB')
xlabel('time,day')
title('Fig.4 Accumulative Water Production vs time')
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')


figure(5)
hold on
plot(T,QGC1);
plot(T,QGC2);
plot(T,QGC3);
plot(T,QGC4);
plot(T,QGC5);
plot(T,QGC6);
plot(T,QGC7);
plot(T,QGC8);
plot(T,QGC9);
plot(T,QGC10);
hold off
grid on
ylabel('Accumulative Gas Production,SCF')
xlabel('time,day')
title('Fig.5 Accumulative Gas Production vs time')
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')

figure(6)
hold on
plot(T,QOC1);
plot(T,QOC2);
plot(T,QOC3);
plot(T,QOC4);
plot(T,QOC5);
plot(T,QOC6);
plot(T,QOC7);
plot(T,QOC8);
plot(T,QOC9);
plot(T,QOC10);
hold off
grid on
ylabel('Accumulative Oil Production,STB')
xlabel('time,day')
title('Fig.6 Accumulative Oil Production vs time')
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')


figure(7)
plot(T(2:73),oil_inc(2:73));
hold on
plot(T(2:73),water_inc(2:73));
plot(T(2:73),0.5*(oil_inc(2:73)+water_inc(2:73)));
hold off
grid on
ylabel('Incremental MB Check')
xlabel('Time,days')
legend('oil phase','water phase','gas phase')
title('Fig.7 Incremental MB Check vs time')

figure(8)
plot(T(2:73),oil_com(2:73));
hold on
plot(T(2:73),water_com(2:73));
plot(T(2:73),0.5*(oil_com(2:73)+water_com(2:73)));
hold off
grid on
ylabel('Cumulative MB Check')
xlabel('Time,days')
legend('oil phase','water phase','gas phase')
title('Fig.8 Cumulative MB Check vs time')
% 
figure(9)
plot(T(2:73),Ro_ab);
hold on
plot(T(2:73),Rw_ab);
plot(T(2:73),Rg_ab);
hold off
grid on
ylabel('Residual')
xlabel('Time,days')
legend('oil phase','water phase','gas phase')
title('Fig.9 Residual vs time')

figure(10)
plot(T,PV1);
hold on
plot(T,PV2);
plot(T,PV3);
plot(T,PV4);
plot(T,PV5);
plot(T,PV6);
plot(T,PV7);
plot(T,PV8);
plot(T,PV9);
plot(T,PV10);
hold off
grid on 
legend('well 1','well 2','well 3','well 4','well 5','well 6','well 7','well 8','well 9','well 10')
ylabel('Well Blocks pressure,psia')
xlabel('time,days')
title('Fig. 10 Well Blocks presure vs time')
