%% -------------------------------------------------------------------------
% Rendezvous01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Rendezvous im Weltall
% Berechnung der Dynamik eines missglückten Weltall Rendezvous 
% von Astronaut außerhalb
% der Raumstation mittels Hill-Gleichung und Clohessy-Wiltshire-Gleichung
% -------------------------------------------------------------------------

%%
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Numerische Lösung
% Parameter
h  = 421;               % Höhe ISS in km
ME = 5.98e24;           % Masse der Erde in kg
RE = 6.378e6;           % Erdradiuse in m
G  = 6.671e-11;         % G in (m^3 / kg /s^2)

x0 = 100;               % Anfangsabstand x in m
y0 = 100;               % Anfangsabstand y in m
z0 = 00;                % Anfangsabstand z in m

v0 = 0.0;               % Anfangsgeschwindigkeit dx in m/s
dx0 = -v0/sqrt(2);      % Anfangsgeschwindigkeit dx in m/s
dy0 = -v0/sqrt(2);      % dto
dz0 = 0;                % dto

R0 = h*1000 + RE;       % Umlaufradius der Raumstation (Target)
w0 = sqrt(G*ME/R0^3);
T0  = 2*pi/w0;
T0h = T0/60;            % Umlaufzeit der Raumstation (Target)

%% 
% Anfangswerte
Y0 = [x0 dx0 y0 dy0 z0 dz0]; % Anfangsvektor

tmax  = 10800;               % Simulationszeit in s
tspan = linspace(0,tmax,1200);
t= tspan;

%% Numerische Lösung volle Lagrangegleichung Linear
opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6);

% Numerische Lösung Hill-Gleichung
[t,Y]=ode45(@dgl_HillGleichung,tspan,Y0,opt,w0); 

x  = Y(:,1);
dx = Y(:,2);
y  = Y(:,3);
dy = Y(:,4);
z  = Y(:,5);
dz = Y(:,6);
kend = 6;
for k=1:kend
    v0 = 0.5-(k-1)*0.1;     % Anfangsgeschwindigkeit v in m/s
    lgdstr(k,:)=num2str(v0,' v_0 = %3.1f m/s');
    dx0 = -v0/sqrt(2);      % Anfangsgeschwindigkeit dx/dt in m/s
    dy0 = -v0/sqrt(2);      % Anfangsgeschwindigkeit dy/dt in m/s
    dz0 = 0;                % Anfangsgeschwindigkeit dz/dt in m/s

    % Lösung Clohessy-Wiltshire-Gleichung
    xCW(k,:)  = (6*y0+4*dx0/w0)*sin(w0*t)+2*dy0*cos(w0*t)/w0-2*dy0/w0 ...
            -(3*dx0+6*w0*y0)*t + x0;
    yCW(k,:)  = -(3*y0+2*dx0/w0)*cos(w0*t)+dy0*sin(w0*t)/w0+4*y0+2*dx0/w0;
    zCW(k,:)  = z0*cos(w0*t)+dz0*sin(w0*t)/w0;
end
% Graphische Ausgabe
figure(1)
hold on
for k=1:kend
    hp(k)=plot(xCW(k,:), yCW(k,:),'color', Colors(2*k-1,:), 'LineWidth',2);
end
ylabel ('\it y \rm in m')
xlabel ('\it x \rm in m')
plot(0,0,'s','MarkerSize',10,'color',Colors(3,:),'LineWidth',3);
plot(x0,y0,'d','MarkerSize',5,'color',Colors(10,:),'LineWidth',3);
grid on
legend(hp,lgdstr,'location','southwest','Numcolumns',2)
legend box off
h=title('CW-Gleichung Bahnebene Astronaut mit v_0 (direkte Linie)');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%% Funktionen
% Hill-Gleichungen
function dY = dgl_HillGleichung(t,Y,w0)
    dY = [   Y(2)
            -2*w0*Y(4)
             Y(4)
            +2*w0*Y(2) + 3*w0^2*Y(3) 
             Y(6)
            -w0^2*Y(5)
         ];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
