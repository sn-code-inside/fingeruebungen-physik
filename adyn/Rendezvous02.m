%% -------------------------------------------------------------------------
% Rendezvous01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Rendezvous im Weltall
% Berechnung der Dynamik eines Weltall-Rendezvous eines Astronauten, der 
% sich außerhalb der Raumstation befindet.
% Hill-Gleichung und Clohessy-Wiltshire-Gleichung
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
RE = 6.378e6;           % Erdradius in m
G  = 6.671e-11;         % G in (m^3 / kg /s^2)

R0 = h*1000 + RE;       % Umlaufradius der Raumstation (Target)
w0 = sqrt(G*ME/R0^3);
T0  = 2*pi/w0;    
T0h = T0/60;            % Umlaufperiode der Raumstation (Target)

%% 

tmax  = 480;           % Simulationszeit  8 min in s
tspan = linspace(0,tmax,120);
t= tspan;

opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6);

%% Numerische Lösung Hill-Gleichung und CW-Gleichungen

% Anfangswerte (Position und Geschwinddigkeit Astronaut(JÃ¤gers)) in m,m/s
XA =250; YA = XA; v0 = 1; 
kend = 5;
for k=1:kend
    x0 = XA-(k-1)*50-25;  
    y0 = x0;                % Anfangsabstand y in m
    lgdstr1(k,:)=num2str(x0,' x_0=y_0 = %3u m');
    z0 = 00;                % Anfangsabstand z in m
    % Astronaut zielt direkt auf Raumstation
    dx0 = -v0/sqrt(2);      % Anfangsgeschwindigkeit dx in m/s
    dy0 = -v0/sqrt(2);      % dto
    dz0 = 0;                % dto
    
    Y0 = [x0 dx0 y0 dy0 z0 dz0]; % Anfangsvektor   
    [t,Y]=ode45(@dgl_HillGleichung,tspan,Y0,opt,w0); 
    
    x(k,:)   = Y(:,1);
    dx(k,:)  = Y(:,2);
    y(k,:)   = Y(:,3);
    dy(k,:)  = Y(:,4);
    z(k,:)   = Y(:,5);
    dz(k,:)  = Y(:,6);
    
    % Lösung Clohessy-Wiltshire-Gleichung  
    tflight = sqrt(x0^2+y0^2)/v0;
    sinw0 = sin(w0*tflight); cosw0 = cos(w0*tflight);
    C       = 3*w0*tflight*sinw0-8*(1-cosw0);
    dx0      = (w0*x0*sinw0-w0*y0*(6*w0*tflight*sinw0-14*(1-cosw0)))/C;
    dy0      = (-2*w0*x0*(1-cosw0)+w0*y0*(4*sinw0-3*w0*tflight*cosw0))/C;
    v0       = sqrt(dx0^2+dy0^2);
    alpha    = atand(dy0/dx0);
    lgdstr2(k,:)=num2str(alpha,' alpha = %3.1f Â°');

    xCW(k,:)  = (6*y0+4*dx0/w0)*sin(w0*t)+2*dy0*cos(w0*t)/w0-2*dy0/w0 ...
            -(3*dx0+6*w0*y0)*t + x0;
    yCW(k,:)  = -(3*y0+2*dx0/w0)*cos(w0*t)+dy0*sin(w0*t)/w0+4*y0+2*dx0/w0;
    zCW(k,:)  = z0*cos(w0*t)+dz0*sin(w0*t)/w0;
end

% Graphische Ausgabe
figure(1)
subplot(121)
hold on
for k=1:kend 
    hp(k)=plot(x(k,:), y(k,:),'color', Colors(2*k-1,:), 'LineWidth',2);
    plot(x(k,1),y(k,1),'d','MarkerSize',5,'color',...
         Colors(2*k-1,:),'LineWidth',3);
end
axis equal
axis square

ylabel ('\it y \rm in m')
xlabel ('\it x \rm in m')
axis ([0 XA -100 YA])           %windows size
grid on
plot(0,0,'s','MarkerSize',10,'color',Colors(3,:),'LineWidth',3);
grid on
legend(hp,lgdstr1,'location','southeast','Numcolumns',1)
legend box off
h=title('Hill-Gleichung direkter Anflug');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(122)
hold on
for k=1:kend 
    hp(k)=plot(xCW(k,:), yCW(k,:),'color', Colors(2*k-1,:), 'LineWidth',2);
    plot(xCW(k,1),yCW(k,1),'d','MarkerSize',5,'color',...
         Colors(2*k-1,:),'LineWidth',3);
end
axis equal
axis square
axis ([0 XA -100 YA])           %windows size
ylabel ('\it y \rm in m')
xlabel ('\it x \rm in m')
grid on 
plot(0,0,'s','MarkerSize',10,'color',Colors(3,:),'LineWidth',3);
grid on
legend(hp,lgdstr2,'location','southeast','Numcolumns',1)
legend box off
h=title('CW-Gleichungen Korrektur im Anflug');
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
