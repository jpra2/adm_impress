% foram plotados para apresentção de CILAMCE 
% problema GAO e Wu 2010, sobre malha triangular distorcida distorcida
figure(1)
MPFA_H=[log2(1/0.125),  log2(0.0027);
        log2(1/0.0625), log2(7.4151e-4);
        log2(1/0.03125), log2(1.9091e-4);
        log2(1/0.015625), log2(4.7358e-5)];
plot(MPFA_H(:,1),MPFA_H(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.125),  log2(0.00102119);
        log2(1/0.0625), log2(2.9987e-4);
        log2(1/0.03125), log2(7.1262e-5);
        log2(1/0.015625), log2(1.8449e-5)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.125),  log2(0.001427);
        log2(1/0.0625), log2(4.1076e-4);
        log2(1/0.03125), log2(1.0659e-4);
        log2(1/0.015625), log2(2.5092e-5)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')
grid
legend('MPFA-H','MPFA-O','MPFA-FPS')
ylabel('log2(Ep(h))')
xlabel('log2(1/h)')

figure(2)

MPFA_H=[log2(1/0.125),  log2(0.0191);
        log2(1/0.0625), log2(0.01);
        log2(1/0.03125), log2(0.0042);
        log2(1/0.015625), log2(0.0021)];
plot(MPFA_H(:,1),MPFA_H(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.125),  log2(0.011068);
        log2(1/0.0625), log2(0.00666);
        log2(1/0.03125), log2(0.002515);
        log2(1/0.015625), log2(0.00131)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.125),  log2(0.0033976);
        log2(1/0.0625), log2(0.001636);
        log2(1/0.03125), log2(4.79384e-4);
        log2(1/0.015625), log2(2.5906e-4)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')
grid
legend('MPFA-H','MPFA-O','MPFA-FPS')
ylabel('log2(Ef(h))')
xlabel('log2(1/h)')
%% foram plotados para apresentção de SPE 
% problema GAO e Wu 2010, sobre malha Kershaw quadrilateral e Triangular
figure(3)
% triangular
MPFA_HD=[log2(1/0.08333),  log2(0.0064);
        log2(1/0.041666), log2(0.0033);
        log2(1/0.0208333), log2(0.0013);
        log2(1/0.01041666), log2(3.7355e-4);
        log2(1/0.00520833), log2(9.9038e-5)];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.0094631);
        log2(1/0.041666), log2(0.004161);
        log2(1/0.0208333), log2(0.0023682);
        log2(1/0.01041666), log2(3.6649e-4);
        log2(1/0.00520833), log2(1.0475e-4)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.034464);
        log2(1/0.041666), log2(0.02795);
        log2(1/0.0208333), log2(0.033498);
        log2(1/0.01041666), log2(0.0014453);
        log2(1/0.00520833), log2(2.3032e-4)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')

hold on

reference=[log2(1/0.08333),  log2(3.2185e-4);
        log2(1/0.041666), log2(8.0565e-5);
        log2(1/0.0208333), log2(2.0144e-5);
        log2(1/0.01041666), log2(5.0359e-6);
        log2(1/0.00520833), log2(1.2588e-6)];

plot(reference(:,1),reference(:,2),'k-')

grid
legend('MPFA-HD','MPFA-O','MPFA-FPS','Reference')
ylabel('log2(Ep(h))')
xlabel('log2(1/h)')

figure(4)
% plotando velocidade
MPFA_HD=[log2(1/0.08333),  log2(0.0903);
        log2(1/0.041666), log2(0.0406);
        log2(1/0.0208333), log2(0.0161);
        log2(1/0.01041666), log2(0.0056);
        log2(1/0.00520833), log2(0.002)];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.1239);
        log2(1/0.041666), log2(0.05488);
        log2(1/0.0208333), log2(0.03125);
        log2(1/0.01041666), log2(0.005336);
        log2(1/0.00520833), log2(0.003205)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.49957);
        log2(1/0.041666), log2(0.715464);
        log2(1/0.0208333), log2(1.8314858);
        log2(1/0.01041666), log2(0.1537);
        log2(1/0.00520833), log2(0.04448)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')
hold on
reference=[log2(1/0.08333),  log2(0.001613);
        log2(1/0.041666), log2(4.3743e-4);
        log2(1/0.0208333), log2(1.1664e-4);
        log2(1/0.01041666), log2(3.0772e-5);
        log2(1/0.00520833), log2(8.0612e-6)];

plot(reference(:,1),reference(:,2),'k-')
grid
legend('MPFA-HD','MPFA-O','MPFA-FPS','Reference')
ylabel('log2(Ef(h))')
xlabel('log2(1/h)')
%%
% plotando Malha Kershaw Quadrilateral
figure(5)
% triangular
MPFA_HD=[log2(1/0.08333),  log2(0.0037);
        log2(1/0.041666), log2(0.0017);
        log2(1/0.0208333), log2(6.2554e-4);
        log2(1/0.01041666), log2(1.8087e-4);
        log2(1/0.00520833), log2(4.7985e-5)];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.00108);
        log2(1/0.041666), log2(2.7616e-4);
        log2(1/0.0208333), log2(6.9651e-5);
        log2(1/0.01041666), log2(1.7474e-5);
        log2(1/0.00520833), log2(4.3757e-6)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.0026299);
        log2(1/0.041666), log2(6.8402e-4);
        log2(1/0.0208333), log2(1.7515e-4);
        log2(1/0.01041666), log2(4.4107e-5);
        log2(1/0.00520833), log2(1.1039e-5)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')
hold on

Reference=[log2(1/0.08333),  log2(9.7417e-4);
        log2(1/0.041666), log2(2.2440e-4);
        log2(1/0.0208333), log2(6.1049e-5);
        log2(1/0.01041666), log2(1.5265e-5);
        log2(1/0.00520833), log2(3.8166e-6)];

plot(Reference(:,1),Reference(:,2),'k-')

grid
legend('MPFA-HD','MPFA-O','MPFA-FPS','Reference')
ylabel('log2(Ep(h))')
xlabel('log2(1/h)')

figure(6)
% plotando velocidade
MPFA_HD=[log2(1/0.08333),  log2(0.0838);
        log2(1/0.041666), log2(0.0321);
        log2(1/0.0208333), log2(0.0113);
        log2(1/0.01041666), log2(0.004);
        log2(1/0.00520833), log2(0.0015)];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.010049);
        log2(1/0.041666), log2(0.00284);
        log2(1/0.0208333), log2(7.7536e-4);
        log2(1/0.01041666), log2(2.0768e-4);
        log2(1/0.00520833), log2(5.5077e-5)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.02255);
        log2(1/0.041666),   log2(0.007106);
        log2(1/0.0208333),  log2(0.002074);
        log2(1/0.01041666), log2(5.8851e-4);
        log2(1/0.00520833), log2(1.6984e-4)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')

hold on 
Reference=[log2(1/0.08333),  log2(0.003082);
        log2(1/0.041666), log2(8.5449e-4);
        log2(1/0.0208333), log2(2.3082e-4);
        log2(1/0.01041666), log2(6.1385e-5);
        log2(1/0.00520833), log2(1.6164e-5)];

plot(Reference(:,1),Reference(:,2),'k-')
grid
legend('MPFA-HD','MPFA-O','MPFA-FPS','Reference')
ylabel('log2(Ef(h))')
xlabel('log2(1/h)')
%% MALHA KERSHAW TRIANGULAR DISTORCIDA
figure(7)
% triangular
MPFA_HD=[log2(1/0.08333),  log2(0.0622);
        log2(1/0.041666), log2(0.0332);
        log2(1/0.0208333), log2(0.0141);
        log2(1/0.01041666), log2(0.0048);
        log2(1/0.00520833), log2(0.0014)];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.98857);
        log2(1/0.041666), log2(0.26642);
        log2(1/0.0208333), log2(0.3512);
        log2(1/0.01041666), log2(0.2487);
        log2(1/0.00520833), log2(0.2555)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.255343201285273);
        log2(1/0.041666), log2(0.253083036249697);
        log2(1/0.0208333), log2(0.254683150707859);
        log2(1/0.01041666), log2(0.255035389628951);
        log2(1/0.00520833), log2(0.259621134567321)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')
hold on
grid
legend('MPFA-HD','MPFA-O','MPFA-FPS')
ylabel('log2(Ep(h))')
xlabel('log2(1/h)')

figure(8)
% plotando velocidade
MPFA_HD=[log2(1/0.08333),  log2(0.1543);
        log2(1/0.041666), log2(0.0757);
        log2(1/0.0208333), log2(0.030);
        log2(1/0.01041666), log2(0.010);
        log2(1/0.00520833), log2(0.0043)];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.6272);
        log2(1/0.041666), log2(0.3159);
        log2(1/0.0208333), log2(0.4489);
        log2(1/0.01041666), log2(0.08838);
        log2(1/0.00520833), log2(0.09123)];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.208259308268331);
        log2(1/0.041666),   log2(0.180889072748497);
        log2(1/0.0208333),  log2(0.156394996448212);
        log2(1/0.01041666), log2(0.109267216116440);
        log2(1/0.00520833), log2(0.094521342123321)];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')

grid
legend('MPFA-HD','MPFA-O','MPFA-FPS')
ylabel('log2(Ef(h))')
xlabel('log2(1/h)')
%% MALHA QUADRILATERAL DISTORCIDA
figure(9)
MPFA_HD=[log2(1/0.08333),  log2(0.0091);
        log2(1/0.041666), log2(0.0060);
        log2(1/0.0208333), log2(5.2789e-04);
        log2(1/0.01041666), log2();
        log2(1/0.00520833), log2()];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.253704384303313);
        log2(1/0.041666), log2(0.254675428020385);
        log2(1/0.0208333), log2(0.254938287510401);
        log2(1/0.01041666), log2(0.255014102606362);
        log2(1/0.00520833), log2()];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.255175708867753);
        log2(1/0.041666), log2(0.255202672328997);
        log2(1/0.0208333), log2(0.255028998050173);
        log2(1/0.01041666), log2(0.255031511840828);
        log2(1/0.00520833), log2()];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')
hold on
grid
legend('MPFA-HD','MPFA-O','MPFA-FPS')
ylabel('log2(Ep(h))')
xlabel('log2(1/h)')

figure(10)
% plotando velocidade
MPFA_HD=[log2(1/0.08333),  log2(0.0303);
        log2(1/0.041666), log2(0.0118);
        log2(1/0.0208333), log2(0.0036);
        log2(1/0.01041666), log2();
        log2(1/0.00520833), log2()];
plot(MPFA_HD(:,1),MPFA_HD(:,2),'b-o')
hold on

MPFA_O=[log2(1/0.08333),  log2(0.025595273889032);
        log2(1/0.041666), log2(0.017885342347609);
        log2(1/0.0208333), log2(0.017222676232888);
        log2(1/0.01041666), log2(0.017173119739437);
        log2(1/0.00520833), log2()];

plot(MPFA_O(:,1),MPFA_O(:,2),'r-<')
hold on 
MPFA_FPS=[log2(1/0.08333),  log2(0.269983132766127);
        log2(1/0.041666),   log2(0.020084157807103);
        log2(1/0.0208333),  log2(0.017430126093869);
        log2(1/0.01041666), log2(0.017228804407073);
        log2(1/0.00520833), log2()];

plot(MPFA_FPS(:,1),MPFA_FPS(:,2),'g-s')

grid
legend('MPFA-HD','MPFA-O','MPFA-FPS')
ylabel('log2(Ef(h))')
xlabel('log2(1/h)')
