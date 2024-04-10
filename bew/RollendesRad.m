% -------------------------------------------------------------------------
% RollendesRad.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Lösung der Bewegungsgleichungen des rollenden Rades.
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = {"-";"--";":"};

%% Variablen:
g = 9.81; % Fallbeschleunigung auf der Erde
R = 0.15; % Radius der Scheibe

% Integrationszeit
tmax = 3.; % Zahl der Sekunden, die betrachtet werden sollen
% Startwerte fuer die Winkel
phi   = 0.;
theta = 3*pi/8;
psi   = 0.;
% Startwerte fuer die Winkelgeschwindigkeiten (rad/s)
omega_phi   = 5.;
omega_theta = 0.;
omega_psi   = 6.;
% Startwerte fuer den Auflagepunkt
x0 = 0.;
y0 = 0.;

%% Numerische Berechnung
tv      = linspace(0.,tmax,10000);
Y0      = [phi;theta;psi;omega_phi;omega_theta;omega_psi;x0;y0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 

[t1,Y1] = ode45(@DGL_Scheibe,tv,Y0,options,g,R);

[t2,Y2] = ode45(@DGL_Reifen,tv,Y0,options,g,R);

%% Graphische Ausgabe
figure()
subplot(2,1,2)
plot(t1,Y1(:,1)/pi,'Linewidth',2,'Color',Colors(2,:),'LineStyle',Style{1});
hold on
plot(t2,Y2(:,1)/pi,'Linewidth',2,'Color',Colors(2,:),'LineStyle',Style{3});
ylabel(['$\varphi/\pi$, $\psi/\pi$'],'Interpreter','Latex', ...
    'FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on;
plot(t1,Y1(:,3)/pi,'Linewidth',2,'Color',Colors(4,:),'LineStyle',Style{1});
plot(t2,Y2(:,3)/pi,'Linewidth',2,'Color',Colors(4,:),'LineStyle',Style{3});
grid on;
legend('$\varphi$ Kreisscheibe','$\varphi$ Reifen','$\psi$ Kreisscheibe',...
    '$\psi$ Reifen','Interpreter','Latex','FontSize',14, ...
    'Location','northwest');
legend box off
set(gca,'FontSize',16);  
hold off;

subplot(2,1,1)
plot(t1,Y1(:,2)/pi,'Linewidth',2,'Color',Colors(2,:),'LineStyle',Style{1});
hold on
plot(t2,Y2(:,2)/pi,'Linewidth',2,'Color',Colors(7,:),'LineStyle',Style{3});
ylabel('\vartheta/\pi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
axis ([0,tmax,0,1.25*max(Y1(:,2)/pi)]);
legend('Kreisscheibe','Reifen','FontSize',14, ...
    'Location','northwest');
legend box off
grid on;
set(gca,'FontSize',16); 

mx = max(Y2(:,7));
my = max(Y2(:,8));
figure()
plot(Y1(:,7),Y1(:,8),'Linewidth',2,'Color',Colors(7,:),'LineStyle',Style{1});
hold on
axis([-1, 1, -0.5, 1.1*my]);
axis equal
plot(Y2(:,7),Y2(:,8),'Linewidth',2,'Color',Colors(9,:),'LineStyle',Style{3});
ylabel('{\it y_A} \rm in m ','FontSize',14);
xlabel('{\it x_A} \rm in m ','FontSize',14);
legend('Kreisscheibe','Reifen','FontSize',14, ...
    'Location','northwest');
legend box off
grid on;
set(gca,'FontSize',16);


% Simulation
figure()
hold on
box on
grid on;
axis equal;
axis([-1, 1, -0.5, 1.1*my]);
set(gca,'FontSize',16);
ylabel('{\it y_A} \rm in m ','FontSize',14);
xlabel('{\it x_A} \rm in m ','FontSize',14);

for ki=1:20:length(t1)
     hold on
     strtime = strjoin({' t = ', num2str(t1(ki), '% 6.2f s')});
     hs(1)= text(-0.9, my ,strtime);
     hs(2)=plot(Y1(ki,7),Y1(ki,8), 'o', 'MarkerFaceColor', Colors(2,:), ...
       'MarkerEdgeColor', Colors(2,:), 'MarkerSize', 2);
     hs(3)=plot(Y2(ki,7),Y2(ki,8), 'o', 'MarkerFaceColor', Colors(3,:), ...
       'MarkerEdgeColor', Colors(3,:), 'MarkerSize', 2);
     pause(.01)
     hs(1).Visible = 'off';
end
text(-0.9, my ,strtime);
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Differentialgleichngen
function dYdt = DGL_Scheibe(t,Y,g,R)
    % Y(1) = phi
    % Y(2) = theta
    % Y(3) = psi
    % Y(4) = d(phi)/dt
    % Y(5) = d(theta)/dt
    % Y(6) = d(psi)/dt
    % Y(7) = x
    % Y(8) = y
    dYdt = [Y(4)
            Y(5)
            Y(6)
            2*Y(5)*Y(6)/sin(Y(2))
            -sin(Y(2))*cos(Y(2))*Y(4)^2-6*sin(Y(2))*Y(4)*Y(6)/5.0-4*g/R ...
              *cos(Y(2))/5.
            5*sin(Y(2))*Y(4)*Y(5)/3.-2*cos(Y(2))*Y(5)*Y(6)/sin(Y(2))
            -R*cos(Y(1))*Y(6)
            -R*sin(Y(1))*Y(6)];
end

function dYdt = DGL_Reifen(t,Y,g,R)
    % Y(1) = phi
    % Y(2) = theta
    % Y(3) = psi
    % Y(4) = d(phi)/dt
    % Y(5) = d(theta)/dt
    % Y(6) = d(psi)/dt
    % Y(7) = x
    % Y(8) = y
    dYdt = [Y(4)
            Y(5)
            Y(6)
            2*Y(5)*Y(6)/sin(Y(2))
            -sin(Y(2))*cos(Y(2))*Y(4)^2-4*sin(Y(2))*Y(4)*Y(6)/2.0-2*g/R ...
              *cos(Y(2))/3.
            3*sin(Y(2))*Y(4)*Y(5)/2.-2*cos(Y(2))*Y(5)*Y(6)/sin(Y(2))
            -R*cos(Y(1))*Y(6)
            -R*sin(Y(1))*Y(6)];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
