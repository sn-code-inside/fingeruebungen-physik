% -------------------------------------------------------------------------
% SchwererKreiselPraez02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Beispielrechnungen und Simulation zur Präzession des 
% schweren symmetrischen Kreisels.
% 
% 
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", "--", "-", ":"];

%% Variablen: Traegheitsmomente, Masse
%
J3     = 1.0;
J12    = 2*J3/3;
omega3 = 10.;
m      = 1.;
g      = 9.81;
s      = 0.2;

% Integrationszeit in s
tmax    = 5.;
% Komponenten der Rotationsachse im koerperfesten System zu t=0
phi0   =  0.;
theta0 =  pi/3.;
psi0   =  0.;
% Erhaltungsgroessen
pphi    = 0.5;
omega3  = 1.;
dtheta0 = 2.25;

%% Potential für theta
%
% Erstelle Graph des Potentials
% [xx,yy]=fplot(@(x)EffektivesPotential(x,pphi,omega3,J12,J3,m,g,s)/...
%               J3/omega3^2,[0.05 pi-0.05]);
theta     = linspace(0.01,pi-0.1,500);
Pot_theta = EffektivesPotential(theta,pphi,omega3,J12,J3,m,g,s)/J3/omega3^2;
figure()
plot(rad2deg(theta),Pot_theta,'Linewidth',1,'Color',Colors(2,:),...
    'LineStyle',Style{1})
axis([0 180 0 5.]);
hold on
% Berechne und zeichne verschiedene Gesamtenergien ein
for k=1:3
    dtheta0k(k) = dtheta0*((k-1)/2)+dtheta0/4;
    Eges(k) = (EffektivesPotential(theta0,pphi,omega3,J12,J3,m,g,s)...
       +0.5*J12*dtheta0k(k)^2)/J3/omega3^2;
    for km=1:length(Pot_theta)-1
       if (Pot_theta(km)-Eges(k)) > 0 && (Pot_theta(km+1)-Eges(k)) < 0 
          m1(k) = km+1;
       else
         if (Pot_theta(km)-Eges(k)) < 0 && (Pot_theta(km+1)-Eges(k)) > 0
           m2(k) = km;
         end
       end
    end
    plot([rad2deg(theta(m1(k))),rad2deg(theta(m2(k)))],[Eges(k),Eges(k)],...
         'Linewidth',1,'Color',Colors(2*k+1,:),'LineStyle',Style{k});
end
ylabel('{\it V}_{\rm eff}(\vartheta)/{\it J}_3\omega_3^2','FontSize',14);
xlabel('\vartheta in °','FontSize',14);
set(gca,'FontSize',14);
grid on

% Suche das Minimum fuer die reine Praezessionsbewegung und uberpruefe es
theta0 = fminsearch(@(x)EffektivesPotential(x,pphi,omega3,J12,J3,m,g,s),1.);
fprintf('\n theta_0 = ');
theta0d = rad2deg(theta0);
fprintf(num2str(theta0d,'%4.1f')');
fprintf('°\n');
dtheta0=0;  % Anfangsgeschwindigkeit 0
Eges(4) = (EffektivesPotential(theta0,pphi,omega3,J12,J3,m,g,s)...
       +0.5*J12*dtheta0^2)/J3/omega3^2;

% Zeichne Gesamtenergie für theta0 ein
plot([theta0d-15,theta0d+15],[Eges(4),Eges(4)],'Linewidth',1,'Color',...
            Colors(9,:),'LineStyle',Style{1})
plot([theta0d,theta0d],[0,Eges(4)],'Linewidth',1,'Color',...
            Colors(9,:),'LineStyle',Style{1})
set(gca,'FontSize',14);
hold off
dtheta0k(4)=0;
stE = ["\it E_1","\it E_2","\it E_2","\it E_0"];
for k = 1:4
   str1(k,:)= strcat(stE(:,k),' \rm = ',(num2str(Eges(k),'% 4.2f')));
   str2(k,:)= strcat('  (d\vartheta/dt)_0 = ', (num2str(dtheta0k(k),'% 4.2f')));
   if k <= 3 
       text(theta0d-45, Eges(k)+0.25,strcat(str1(k,:),str2(k,:)),'FontSize',12, ...
           'Color', Colors(2*k+1,:)); 
   else
       text(20, Eges(k),strcat(str1(k,:),str2(k,:)),'FontSize',12, ...
           'Color', Colors(2*k+1,:)); 
   end
end


%% Zeitentwicklung theta
%
figure();
hold on;
theta0 =  pi/3.;
dtheta0k(4) = 0;
tspan=[0.,tmax];
for k=1:4
    if k==4 
        theta0=deg2rad(theta0d);
    end
    Y0=[theta0;dtheta0k(k);phi0;psi0];
    options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
    [t1,Y1] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);
     plot(t1,rad2deg(Y1(:,1)),'Linewidth',1,'Color',...
        Colors(k*2+1,:),'LineStyle',Style{k});
end
set(gca,'FontSize',14);
ylabel('\vartheta in °','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
legend(stE,'location','southwest','NumColumns',4);
legend box off
grid on;
hold off

%% Zeitentwicklung der Eulerschen Winkel
%% Ohne Nutation

% Energieauswahl
kE = 4; 

tspan=[0.,tmax];
Y0=[theta0;dtheta0;phi0;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);

figure
% theta
plot(t1,rad2deg(Y1(:,1)),'Linewidth',1,'Color',Colors(11,:),'LineStyle',Style{1});
ylabel('\phi,\vartheta,\psi','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
% phi
plot(t1,rad2deg(Y1(:,3)),'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
% psi
plot(t1,rad2deg(Y1(:,4)),'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
legend('\phi','\vartheta','\psi in °','FontSize',14,...
       'Location','northwest');
legend box off
text(1,360,strcat(stE(:,kE),'\rm  = ', num2str(Eges(kE),'%4.2f')),'FontSize', 14);
set(gca,'FontSize',14);
hold off
grid on;
axis([0 tmax 0 400]);

%% Zeitentwicklung der Eulerschen Winkel
%% Mit Nutation

% Energieauswahl
kE = 1; 

tspan=[0.,tmax];
Y0=[pi/3.;dtheta0k(kE);phi0;psi0];
options = odeset('AbsTol',1.e-8,'RelTol',1.e-6); 
[t1,Y1] = ode45(@DGL,tspan,Y0,options,pphi,omega3,J12,J3,m,g,s);

figure()
% theta
plot(t1,rad2deg(Y1(:,1)),'Linewidth',1,'Color',Colors(11,:),'LineStyle',Style{1});
ylabel('\phi,\vartheta,\psi in °','FontSize',14);
xlabel('{\it t} \rm in s ','FontSize',14);
hold on
% phi
plot(t1,rad2deg(Y1(:,3)),'Linewidth',1,'Color',Colors(4,:),'LineStyle',Style{1});
% psi
plot(t1,rad2deg(Y1(:,4)),'Linewidth',1,'Color',Colors(3,:),'LineStyle',Style{1});
legend('\phi','\vartheta','\psi','FontSize',14,...
       'Location','northwest');
legend box off
text(1,360,strcat(stE(:,kE),'\rm = ', num2str(Eges(kE),'%4.2f')),'FontSize', 14);
set(gca,'FontSize',14);
hold off
grid on;
axis([0 tmax 0 400]);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

% Differentialgleichung
function dYdt = DGL(t,Y,pphi,omega3,J12,J3,m,g,s)
    % Y(1) = theta
    % Y(2) = d(theta)/dt
    % Y(3) = phi
    % Y(4) = psi
    dYdt = [Y(2)
            -(pphi-J3*omega3*cos(Y(1)))*J3*omega3/(J12^2*sin(Y(1)))+...
            (pphi-J3*omega3*cos(Y(1)))^2*cos(Y(1))/(J12^2*sin(Y(1))^3)+...
            m*g*s*sin(Y(1))/J12
            (pphi-J3*omega3*cos(Y(1)))/(J12*sin(Y(1))^2)
            omega3-(pphi-J3*omega3*cos(Y(1)))*cos(Y(1))/(J12*sin(Y(1))^2)];
end

function Pot = EffektivesPotential(theta,pphi,omega3,J12,J3,m,g,s)
     Pot = (pphi - J3*omega3*cos(theta)).^2./(2.*J12*sin(theta).^2) ...
           + 0.5*J3*omega3^2 + m*g*s*cos(theta);
end 
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
