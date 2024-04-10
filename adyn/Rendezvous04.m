%% -------------------------------------------------------------------------
% Rendezvous03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Rendezvous im Weltall
% Berechnung der Kinematik eines Weltall-Rendezvous eines Raumtransporters
% und der ISS (Lösung über die Kepler-Gleichung bzw. CW-Gleichungen) und
% Anflug mit direkter Sichtlinie.
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
h  = 430;               % Höhe ISS in km
ME = 5.98e24;           % Masse der Erde in kg
RE = 6.378e6;           % Erdradius in m
G  = 6.671e-11;         % G in (m^3 / kg /s^2)

R0 = h*1000 + RE;       % Umlaufradius der ISS (Target)
w0 = sqrt(G*ME/R0^3);
T0  = 2*pi/w0;    
T0min = T0/60;          % Umlaufperiode der Raumstation (Target) in min

ecc = 0.035;            % Elliptizität der Raumtransporter (Jäger) Bahn

%% 

tmax  = T0;           % Simulationszeit in s
tspan = linspace(0,tmax,121);
t= tspan;

%% Numerische Lösung Kepler-Gleichung

% Anfangswerte (Position und Geschwindigkeit ISS(Target)) in m,m/s
REA = RE /1e6; R0A = R0/1e6;  % Normierung auf 10^6 m
XA = 0/1e6; YA = R0/1e6; 
VA = w0*R0A;   % in 10^6 m/s
% Berechnung Bewegung ISS  
% (Achtung Anfangsbedingungen gewählt bei (R0A,0))!!!
Y = R0A*cos(w0*t);
X = R0A*sin(w0*t);
XE = REA*cos(w0*t);  % Erde
YE = REA*sin(w0*t);

% Anfangswerte (Position und Geschwinddigkeit Raumschiff(Jägers)) in m,m/s
Rmin = R0A*(1-ecc);
Rmax = R0A*(1+ecc);
GA   = G/1e18;  % Normierung auf 10^6 m
xA = 0/1e6; yA = Rmin; 
vA = sqrt(ME*GA*(2/Rmin-1/R0A));  % in 10^6 m/s

%% Berechnung Bewegung Raumschiff(Jägers)) Kepler-Gleichung
% Mittlere Anomalie
M = w0*t;
% Berechnung der exzentrischen Anomalie
Ecc=EAnomalie(M,ecc);
% Berechnung der Bahndaten aus Kepler-Gleichung
fac = sqrt((1-ecc)*(1+ecc));
cosE = cos(Ecc);
sinE = sin(Ecc);
dis = R0A.*(1-ecc*cosE);
% Koordinaten Target 
% (Achtung Anfangsbedingungen gewählt bei Perigäum)!!!
y   = R0A.*(cosE-ecc);
x   = R0A*fac.*sinE;

% Relativkoordinaten und Drehung des KOS der ISS
xR = x-X;
yR = y-Y;
zR = 0*yR;

for k=1:length(t)
    rTJ  = [xR(k) yR(k) zR(k)]';
    % Drehung um Uhrzeigersinn !
    rTJ2 = mtimes(R_z(-w0*t(k)),rTJ);
    xR(k) = rTJ2(1);
    yR(k) = rTJ2(2);
end



%% Numerische Lösung CW-Gleichung für Relativkoordinaten

% AB
y0 = (y(1)-Y(1));
x0 = 0;
dx0 = -2*y0*w0;
dy0 = 0;
dx01 = -2*y0*sqrt(GA*ME/R0A^3);
ddx0 =dx0-dx01

% CW-Gleichungen
xRnum  = (6*y0+4*dx0/w0)*sin(w0*t)+2*dy0*cos(w0*t)/w0-2*dy0/w0 ...
            -(3*dx0+6*w0*y0)*t + x0;
yRnum  = -(3*y0+2*dx0/w0)*cos(w0*t)+dy0*sin(w0*t)/w0+4*y0+2*dx0/w0;



%% Graphische Ausgabe

figure(1)
subplot(121)
hp(1) = plot(X,Y,'+','MarkerIndices',1:10:length(x),...
    'color',Colors(4,:),'LineWidth',2,'LineStyle',Style(1));
zaehl = 0;
for k=1:10:length(X)-1
    zaehl = zaehl +1;
    hl(k)=LabelPoints(X(k), Y(k), ...
        num2str(zaehl),'N',0.25,'FontSize',14,'Color',Colors(4,:));
end
hold on
hp(2) = plot(x,y,'d','MarkerIndices',1:10:length(x),...
    'color',Colors(2,:),'LineWidth',2,'LineStyle',Style(3));
zaehl = 0;
for k=1:10:length(X)-1
    zaehl = zaehl +1;
    hl(k)=LabelPoints(x(k), y(k), ...
        num2str(zaehl),'S',0.25,'FontSize',14,'Color',Colors(2,:));
end
axis equal
axis square
ylabel ('\it y \rm in 1000 km')
xlabel ('\it x \rm in 1000 km')
axis ([-10 10 -10 10])           %windows size
grid on
hp(3)= plot(XE,YE,'color',Colors(3,:),'LineWidth',2);
legend(hp,'ISS (Target)','Transporter (Jäger)','Erde',...
       'location','south','Numcolumns',2)
legend box off
h=title('Bahn Target und Jäger');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(122)
hp(4) = plot(xR*1e3,yR*1e3,'+','MarkerIndices',1:10:length(x),...
    'color',Colors(2,:),'LineWidth',2,'LineStyle',Style(1));
zaehl = 0;
hold on
for k=1:10:length(X)-1
    zaehl = zaehl +1;
    hl(k)=LabelPoints(xR(k)*1e3, yR(k)*1e3, ...
        num2str(zaehl),'S',0.25,'FontSize',14,'Color',Colors(2,:));
end
hp(5) = plot(xRnum*1e3,yRnum*1e3,'+','MarkerIndices',1:10:length(x),...
    'color',Colors(3,:),'LineWidth',2,'LineStyle',Style(1));
zaehl = 0;
for k=1:10:length(X)-1
    zaehl = zaehl +1;
    hl(k)=LabelPoints(xRnum(k)*1e3, yRnum(k)*1e3, ...
        num2str(zaehl),'N',0.25,'FontSize',14,'Color',Colors(3,:));
end
axis equal
axis square
ylabel ('\it y´ \rm in km')
xlabel ('\it x´ \rm in km')
axis ([-1 1 -1 1]*1.0e3)           %windows size
grid on
hp(6) = plot(0,0,'s','MarkerSize',10,'color',Colors(4,:),'LineWidth',3);
legend(hp(4:6),'Kepler-Gleichung','CW-Gleichungen','ISS',...
       'location','south','Numcolumns',1)
legend box off
h=title('Bahn Jäger von ISS aus');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


%% Numerische Lösung CW-Gleichungen für Direktanflug

% Parameter
h  = 421;               % Höhe ISS in km
ME = 5.98e24;           % Masse der Erde in kg
RE = 6.378e6;           % Erdradius in m
G  = 6.671e-11;         % G in (m^3 / kg /s^2)
R0 = h*1000 + RE;       % Umlaufradius der Raumstation (Target)
w0 = sqrt(G*ME/R0^3);
tmax  = 200;           % Simulationszeit  8 min in s
t     = linspace(0,tmax,201);

% Anfangswerte (Position und Geschwindigkeit Astronaut(Jägers)) in m,m/s
XA =120; v0 = 1; 
kend = 6;
for k=1:kend
    x0 = XA-(k-1)*20;  
    y0 = 0;                % Anfangsabstand y in m
    z0 = 0;                % Anfangsabstand z in m
    % Astronaut zielt direkt auf Raumstation
    dx0 = -v0;             % Anfangsgeschwindigkeit dx in m/s
    dy0 = 0;               % dto
    dz0 = 0;               % dto    
    % Lösung Clohessy-Wiltshire-Gleichung  
    xCW(k,:)  = (6*y0+4*dx0/w0)*sin(w0*t)+2*dy0*cos(w0*t)/w0-2*dy0/w0 ...
            -(3*dx0+6*w0*y0)*t + x0;
    yCW(k,:)  = -(3*y0+2*dx0/w0)*cos(w0*t)+dy0*sin(w0*t)/w0+4*y0+2*dx0/w0;
    zCW(k,:)  = z0*cos(w0*t)+dz0*sin(w0*t)/w0;
end

%% Graphische Ausgabe
figure(2)
hold on
for k=1:kend 
    hp(k)=plot(xCW(k,:), yCW(k,:),'color', Colors(2*k-1,:), 'LineWidth',2);
    plot(xCW(k,1),yCW(k,1),'d','MarkerSize',5,'color',...
         Colors(2*k-1,:),'LineWidth',3);
end
%axis equal
axis ([0 XA -20 5])           %windows size
ylabel ('\it y \rm in m')
xlabel ('\it x \rm in m')
grid on 
plot(0,0,'s','MarkerSize',10,'color',Colors(4,:),...
     'MarkerFaceColor',Colors(4,:),'LineWidth',3);
grid on
h=title('CW-Gleichungen');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%% Funktionen

% Exzentrische Anomalie
function Ecc = EAnomalie(M,e)
    f= 100000;
    E=M;
    while abs(f) > 1e-9
       f=E-e.*sin(E)-M;
       E=E-f./(1-e.*cos(E));
    end
    Ecc=E;
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
