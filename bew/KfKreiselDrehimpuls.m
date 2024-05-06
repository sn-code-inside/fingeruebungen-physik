% -------------------------------------------------------------------------
% KfKreiselDrehimpuls.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Beispielrechnungen und Simulation zur Stabilität des allgemeinen
% kräftefreien Kreisels.
%
% Drehimpuls des allgemeinen kräftefreien Kreisels im körperfesten
% System.
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
tmax = 150;
tspan=[0.,tmax];


%%
% anfaengliche Rotation um zwei Achsen
Y0=[1;0.8;0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);

figure()
subplot(1,2,1)
hp(1) = plot3(Y1(:,1),Y1(:,2),Y1(:,3),'Linewidth',1,'Color',Colors(5,:),...
      'LineStyle',Style(1));
zlabel('\omega_{1} \rm in 1/s ','FontSize',14);
ylabel('\omega_{2} \rm in 1/s ','FontSize',14);
xlabel('\omega_{3} \rm in 1/s ','FontSize',14);
grid on
hold on
for k = 1:2:length(Y1(:,1))
plot3([0; Y1(k,1)], [0; Y1(k,2)], [0; Y1(k,3)],'Color',Colors(5,:),...
      'LineStyle',Style(1));
end
axis([-1 1 -1 1 -1 1]);
xl = [-1 1];
yl = [-1 1];
[X,YB] = meshgrid(xl,yl);
surf(X,YB,zeros(size(X)));
shading flat
alpha 0.1
title('Vektor \omega (t) für Anfangsrotation \omega = (1,0.8,0)', ...
    'FontWeight', 'normal');
legend off
set(gca,'FontSize',14);

subplot(1,2,2)
plot(t1,J1*Y1(:,1),'Linewidth',1,'Color',Colors(2,:),'LineStyle', ...
    Style{1});
%axis([0 tmax -1.1 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
grid on
hold on
plot(t1,J2*Y1(:,2),'Linewidth',1,'Color',Colors(3,:),'LineStyle', ...
    Style{1});
plot(t1,J3*Y1(:,3),'Linewidth',1,'Color',Colors(4,:),'LineStyle', ...
    Style{1});
plot(t1,J1*Y1(:,1).^2+J2*Y1(:,2).^2+J3*Y1(:,3).^2,'Linewidth',2, ...
    'Color',Colors(5,:),'LineStyle', Style{1});
legend('L_{1}','L_{2}','L_{3}','L^2','FontSize',14,...
       'Location','southoutside','Orientation','horizontal');
legend box off
title('Komponenten und Quadrat des Drehimpulses', 'FontWeight', ...
    'normal');
set(gca,'FontSize',14);


%%
% anfaengliche Rotation um alle drei Achsen
Y0=[0.3;0.8;0.4];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,J1,J2,J3);

figure()
subplot(1,2,1)
hp(1) = plot3(Y1(:,1),Y1(:,2),Y1(:,3),'Linewidth',1,'Color',Colors(5,:),...
      'LineStyle',Style(1));
zlabel('\omega_{1} \rm in 1/s ','FontSize',14);
ylabel('\omega_{2} \rm in 1/s ','FontSize',14);
xlabel('\omega_{3} \rm in 1/s ','FontSize',14);
grid on
hold on
for k = 1:2:length(Y1(:,1))
plot3([0; Y1(k,1)], [0; Y1(k,2)], [0; Y1(k,3)],'Color',Colors(5,:),...
      'LineStyle',Style(1));
end
axis([-1 1 -1 1 -1 1]);
xl = [-1 1];
yl = [-1 1];
[X,YB] = meshgrid(xl,yl);
surf(X,YB,zeros(size(X)));
shading flat
alpha 0.1
title('Vektor \omega (t) für Anfangsrotation \omega = (0.3,0.8,0.4)', ...
    'FontWeight', 'normal');
legend off
set(gca,'FontSize',14);

subplot(1,2,2)
plot(t1,J1*Y1(:,1),'Linewidth',1,'Color',Colors(2,:),'LineStyle', ...
    Style{1});
%axis([0 tmax -1.1 1.1]);
ylabel('\omega_{\it i} \rm in 1/s ','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
grid on
hold on
plot(t1,J2*Y1(:,2),'Linewidth',1,'Color',Colors(3,:),'LineStyle', ...
    Style{1});
plot(t1,J3*Y1(:,3),'Linewidth',1,'Color',Colors(4,:),'LineStyle', ...
    Style{1});
plot(t1,J1*Y1(:,1).^2+J2*Y1(:,2).^2+J3*Y1(:,3).^2,'Linewidth',2, ...
    'Color',Colors(5,:),'LineStyle', Style{1});
legend('L_{1}','L_{2}','L_{3}','L^2','FontSize',14,...
       'Location','southoutside','Orientation','horizontal');
legend box off
title('Komponenten und Quadrat des Drehimpulses', 'FontWeight', ...
    'normal');
set(gca,'FontSize',14);
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


