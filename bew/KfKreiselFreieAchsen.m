% -------------------------------------------------------------------------
% KfKreiselFreieAchsen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Stabilität des allgemeinen kräftefreien Kreisels.
% 
% Beispielrechnungen und Simulation zur Stabilität des allgemeinen
% kräftefreien Kreisels.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Variablen: Traegheitsmomente
J1 = 1.5;
J2 = 2;
J3 = 2.5;
% Integrationszeit
tmax = 120;
tspan=[0.,tmax];


%%
% kleine Störung auf Hautträgheitsachse 2, Hauptdrehung um Achse 1
stoerg1k = 0.05;
stoerg1m = 0.25;
stoerg1g = 0.50;
Y0=[1;stoerg1k;0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
Y0=[1;stoerg1m;0];
[t2,Y2] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
Y0=[1;stoerg1g;0];
[t3,Y3] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);


figure()
subplot(1,2,1)
plot(t1,Y1(:,1),'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style{1});
axis([0 tmax -1.1 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
grid on
hold on
plot(t1,Y1(:,2),'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
plot(t1,Y1(:,3),'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
legend('\omega_{1}','\omega_{2}','\omega_{3}','FontSize',14,...
       'Location','best');
legend box off
title(strcat('\Delta \omega = ', num2str(stoerg1k,3), ' 1/s'), 'FontWeight','normal');
set(gca,'FontSize',14);


subplot(1,2,2)
hp(1) = plot3(Y1(:,1),Y1(:,2),Y1(:,3),'Linewidth',1,'Color',Colors(5,:),...
      'LineStyle',Style(1));
zlabel('\omega_{1} \rm in 1/s ','FontSize',14);
ylabel('\omega_{2} \rm in 1/s ','FontSize',14);
xlabel('\omega_{3} \rm in 1/s ','FontSize',14);
grid on
hold on
hp(2) = plot3(Y2(:,1),Y2(:,2),Y2(:,3),'Linewidth',1,'Color',Colors(7,:),...
      'LineStyle',Style(1));
hp(3) = plot3(Y3(:,1),Y3(:,2),Y3(:,3),'Linewidth',1,'Color',Colors(9,:),...
      'LineStyle',Style(1));
for k = 1:2:length(Y1(:,1))
plot3([0; Y1(k,1)], [0; Y1(k,2)], [0; Y1(k,3)],'Color',Colors(5,:),...
      'LineStyle',Style(1));
end
axis([-1+stoerg1g 1+stoerg1g -1 1 -1 1]);
xl = [-1+stoerg1g 1+stoerg1g];
yl = [-1 1];
[X,YB] = meshgrid(xl,yl);
surf(X,YB,zeros(size(X)));
shading flat
alpha 0.1
legend(hp, strcat('\Delta \omega_1 =',num2str(stoerg1k,3)),...
       strcat('\Delta \omega_1 =',num2str(stoerg1m,3)),...
       strcat('\Delta \omega_1 =',num2str(stoerg1g,3)),...
              'FontSize',14,'Location','south');
title(' Vektor \omega (t) für verschiedene Störungen \Delta \omega',...
      'FontWeight','normal');
legend box off
set(gca,'FontSize',14);


%%
% kleine Störung auf Hautträgheitsachse 1, Hauptdrehung um Achse 2
stoerg1k = 0.05;
stoerg1m = 0.25;
stoerg1g = 0.50;
Y0=[stoerg1k;1;0];
[t1,Y1] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
% mitteler Störung auf Hautträgheitsachse 1, Hauptdrehung um Achse 2
Y0=[stoerg1m;1;0];
[t2,Y2] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);
% große Störung auf Hautträgheitsachse 1, Hauptdrehung um Achse 2
Y0=[stoerg1g;1;0];
[t3,Y3] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);

figure
subplot(1,2,1)
plot(t1,Y1(:,1),'Linewidth',1,'Color',Colors(2,:),'LineStyle',Style(1));
axis([0 tmax -1.1 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
grid on
plot(t1,Y1(:,2),'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style(1));
plot(t1,Y1(:,3),'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style(1));
legend('\omega_{1}','\omega_{2}','\omega_{3}','FontSize',14,...
       'Location','best');
title(strcat('\Delta \omega = ', num2str(stoerg1k,3), ' 1/s'), 'FontWeight','normal');
legend box off
set(gca,'FontSize',14);

subplot(1,2,2)
hp(1) = plot3(Y1(:,1),Y1(:,2),Y1(:,3),'Linewidth',1,'Color',Colors(5,:),...
      'LineStyle',Style(1));
zlabel('\omega_{1} \rm in 1/s ','FontSize',14);
ylabel('\omega_{2} \rm in 1/s ','FontSize',14);
xlabel('\omega_{3} \rm in 1/s ','FontSize',14);
grid on
hold on
hp(2) = plot3(Y2(:,1),Y2(:,2),Y2(:,3),'Linewidth',1,'Color',Colors(7,:),...
      'LineStyle',Style(1));
hp(3) = plot3(Y3(:,1),Y3(:,2),Y3(:,3),'Linewidth',1,'Color',Colors(9,:),...
      'LineStyle',Style(1));
for k = 1:2:length(Y1(:,1))
plot3([0; Y1(k,1)], [0; Y1(k,2)], [0; Y1(k,3)],'Color',Colors(5,:),...
      'LineStyle',Style(1));
end
xyzlim=([-1.1 1.1]);
axis([xyzlim(1) xyzlim(2) xyzlim(1) xyzlim(2) xyzlim(1) xyzlim(2)]);
xl = xyzlim;
yl = xyzlim;
[X,YB] = meshgrid(xl,yl);
surf(X,YB,zeros(size(X)))
shading flat
alpha 0.1
legend(hp(1:3), strcat('\Delta \omega_1 =',num2str(stoerg1k,3)),...
       strcat('\Delta \omega_1 =',num2str(stoerg1m,3)),...
       strcat('\Delta \omega_1 =',num2str(stoerg1g,3)),...
              'FontSize',14,'Location','south');
legend box off
set(gca,'FontSize',14);

% print ("AsymmKreiselStabilitaet_instab.pdf","-S652,425","-F:14");
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%%

% Differentialgleichung
function dYdt = DGL(t,Y,J1,J2,J3)
    % Y (i) = omege_i
    dYdt = [(J2-J3)/J1*Y(2)*Y(3)
		    (J3-J1)/J2*Y(3)*Y(1)
		    (J1-J2)/J3*Y(1)*Y(2)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


