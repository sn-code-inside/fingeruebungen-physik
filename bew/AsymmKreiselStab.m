% -------------------------------------------------------------------------
% AsymmKreiselStab.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Stabilität asymmetrischer Kreisel
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = {"-";"--";":"};

% Variablen: Traegheitsmomente
J1 = 1.5;
J2 = 2;
J3 = 2.5;
%Integrationszeit
tmax = 50.;

tspan=[0.,tmax];
Y0=[1;0;0.1];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
Y0=[1;0;0.2];
[t2,Y2] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
Y0=[1;0;0.5];
[t3,Y3] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);

figure
subplot(1,2,1)
plot(t1,Y1(:,1),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
axis([0 tmax -1.1 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
plot(t1,Y1(:,2),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{2});
plot(t1,Y1(:,3),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{3});
legend('\omega_{1}','\omega_{2}','\omega_{3}','FontSize',14,...
       'Location','southwest');
subplot(1,2,2)
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'Linewidth',1,'Color',Colors(1,:),...
      'LineStyle',Style{1});
axis([-1 1 -1 1 -1 1])
zlabel('\omega_{1} \rm in 1/s ','FontSize',14);
ylabel('\omega_{2} \rm in 1/s ','FontSize',14);
xlabel('\omega_{3} \rm in 1/s ','FontSize',14);
grid on
hold on
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'Linewidth',1,'Color',Colors(2,:),...
      'LineStyle',Style{1});
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'Linewidth',1,'Color',Colors(3,:),...
      'LineStyle',Style{1});
axis([-2 2 -2 2 -2 2])
legend('0.1','0.2','0.5','FontSize',14,'Location','northeast');
%print ("AsymmKreiselStabilitaet_stab.pdf","-S652,425","-F:14");
hold off


tspan=[0.,tmax];
Y0=[0.1;1;0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
Y0=[0.2;1;0];
[t2,Y2] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
Y0=[0.5;1;0];
[t3,Y3] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);

figure
subplot(1,2,1)
plot(t1,Y1(:,1),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{1});
axis([0 tmax -1.1 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
plot(t1,Y1(:,2),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{2});
plot(t1,Y1(:,3),'Linewidth',1,'Color',Colors(1,:),'LineStyle',Style{3});
legend('\omega_{1}','\omega_{2}','\omega_{3}','FontSize',14,...
       'Location','southwest');
subplot(1,2,2)
plot3(Y1(:,1),Y1(:,2),Y1(:,3),'Linewidth',1,'Color',Colors(1,:),...
      'LineStyle',Style{1});
axis([-1 1 -1 1 -1 1])
zlabel('\omega_{1} \rm in 1/s ','FontSize',14);
ylabel('\omega_{2} \rm in 1/s ','FontSize',14);
xlabel('\omega_{3} \rm in 1/s ','FontSize',14);
grid on
hold on
plot3(Y2(:,1),Y2(:,2),Y2(:,3),'Linewidth',1,'Color',Colors(2,:),...
      'LineStyle',Style{1});
plot3(Y3(:,1),Y3(:,2),Y3(:,3),'Linewidth',1,'Color',Colors(3,:),...
      'LineStyle',Style{1});
axis([-2 2 -2 2 -2 2])
legend('0.1','0.2','0.5','FontSize',14,'Location','northeast');
%print ("AsymmKreiselStabilitaet_instab.pdf","-S652,425","-F:14");

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%% Funktionen

% Differentialgleichung
function dYdt = DGL(t,Y,J1,J2,J3)
  dYdt = [(J2-J3)/J1*Y(2)*Y(3)
		 (J3-J1)/J2*Y(3)*Y(1)
		 (J1-J2)/J3*Y(1)*Y(2)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


