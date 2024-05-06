% -------------------------------------------------------------------------
% FallenderStab03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik verschiedener fallender Staebe auf Basis der
% Lagrange-Gleichungen  
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerische Lösung

% Parameter
L       =  1.000;           % Länge Stab in  m
m       =  1.000;           % Masse des Stabs in kg
mT      =  5.000;           % Reibungskoeffizient 1 
g       =  9.81;            % g in m/s^2 
theta0  =  deg2rad(1);      % Anfangswinkel 
dtheta0 =  0;               % Anfangsgeschwindigkeit rad/s
tmax    =  2.0;             % Fallzeit in s
tspan   =  linspace(0.0,tmax,200);

fprintf('\n ');
fprintf('\n L   = %8.2f m', L);
fprintf('\n m   = %8.2f kg', m);
fprintf('\n mT   = %8.2f kg', mT);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB=[theta0;dtheta0]; % AB für ode45

%%
% Dynamikberechnung für verschiedene Staebe

% Parameter für DGL 
Para = [3*g/2/L, 33*g/32/L, 54*g/11/L];

Para_Str = string(3);
Para_Str(1,:) = string('homogener Stab');
Para_Str(2,:) = string('Zusatzgewicht am Ende des Stabes');
Para_Str(3,:) = string('Zusatzgewicht symmetrisch um S');

for k=1:3
    P1 = Para(k);
    % Numerische Lösung LGL
    options = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent1);
    [t, Y, TE, YE, IE] = ode45(@dgl_FallingRod1,tspan, AB,options,P1); 
    switch k
    case 1
        t1   = t;
        z1   = rad2deg(Y(:,1));
    case 2
        t2   = t;
        z2   = rad2deg(Y(:,1));
    case 3
        t3   = t;
        z3   = rad2deg(Y(:,1));
    end
end

%% 
% Graphische Ausgabe

tplot = 1.2*max([t1;t2;t3]);
fig=figure();

% Fallwinkel
h=title('Fallende Stäbe');
hold on
lp(1) = plot(t1,L*cosd(z1));
lp(2) = plot(t2,L*cosd(z2));
lp(3) = plot(t3,L*cosd(z3));
for k=1:3 
    set(lp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
end
axis([0 tplot -0.05 1.05]);
grid on
ylabel('Höhe h  in m','FontSize',14)
xlabel('t in s','FontSize',14)
% legend(lp, Para_Str,'location','southwest');
% legend box off
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Functions

function dY = dgl_FallingRod1(t,Y,P1)
    % q- Winkel theta(t), 
    % dq- Winkelgeschwindigkeit dtheta(t)
    dY    = [Y(2);...    
             P1*sin(Y(1))];
end

 
% Ereignisfunktion

function [value,isterminal,direction] = MyEvent1(t,Y,P1)
    % Ereignisfunktion bis Eintreten des Gleitens bei Slip Winkel thetaS bzw.
    % thetaC
        value = (pi/2-Y(1)); % erkenne 90°
        isterminal = 1;      % stop Integration
        direction = 0;       % negative Richtung
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------



