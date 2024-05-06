% -------------------------------------------------------------------------
% Ionantrieb02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet numerisch (ODE45) die Bahnen/Flugzeiten
% für einen Ionenantrieb aus einem LEO in ein GEO 
% im Vergleich zum Hohmannn-Transfer 
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Parameter
GME = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius in km
h1  = 1000;          % Höhe LEO Bahn
r1  = h1+RE;         % geozentrischer Radius Ausgangsbahn
h2  = 35786;         % Höhe GEO Bahn
r2  = RE+h2;
rGEO= RE+h2;

% Raumschiff
m0 = 10000;  % Anfangsmasse in kg
FS = 0.5;    % Schub in N
mdot = 3e-3; % Massenausstoß in kg/s

% Hilfsgrößen
phi = linspace(0,2*pi,361);
xK    = cos(phi);
yK    = sin(phi);

%% Berechnung des Geschwindigkeitszuwachs und der Bahn nach Hohmann
aH0 = 0.5*(r1+r2);            % Transferellipse Halbachse
vP = sqrt(GME*r2/aH0/r1);     % v bei Perigäum
vA = sqrt(GME*r1/aH0/r2);     % v bei Apogäum
tH0= pi*sqrt(aH0^3/GME)/3600; % Flugzeit für Hohmann Ellipse in h
eccH0 = sqrt(1-r1*r2/aH0^2);  % Exzentrizität

tH = linspace(0,tH0*3600,1000);

% Mittlere Anomalie
M = sqrt(GME./aH0.^3).*tH;

% Berechnung der exzentrischen Anomalie
Ecc=EAnom(M,eccH0);
Ecc=Ecc;
cosE=cos(Ecc);
sinE=sin(Ecc);
fac=sqrt((1-eccH0)*(1+eccH0));

% Bahnkoordinaten 
xH     = aH0.*(cosE-eccH0);
yH     = aH0.*fac.*sinE; 
rH     = sqrt(xH.^2+yH.^2);

%% Numerische Berechnung Ionenantrieb Abstand Erde 
tend = 5*tH0*3600;
t  = linspace(0,tend,10000);
x0  = RE+h1;
y0  = 0; 
vx0 = 0; 
vy0 = sqrt(GME/x0);                
AB = [x0;vx0;y0;vy0];          %AB für DGL

P1.m0   = m0;
P1.GME  = GME;
P1.FS   = FS;
P1.mdot = mdot;

% MATLABs Runge-Kutta ode45 Routine 
opts = odeset('AbsTol',1.e-9,'RelTol',1.e-8);
[tS,Y]=ode45(@(t,Y,P1)DGL_IonTriebwerk(t,Y,P1),[0 tend],AB,opts,P1);

xS = Y(:,1);
yS = Y(:,3);
rS = sqrt(xS.^2+yS.^2);
ratioF = FS/m0./(GME./rS.^2);

for k=1:length(rS)
    if rS(k) > r2
        rS(k) = NaN;
        xS(k) = NaN;
        ratioF(k) = NaN;
    end
end
[maxrS, kend] = max(rS);

%%

figure('Name','Ionenantrieb Abstand Erde und Gravitation/Schub über Zeit')
yyaxis left
h(1)=plot(tS/3600, rS/1000,'color',Colors(3,:),'Linewidth',2,'LineStyle',Style(1));
hold on
h(2)=plot(tH/3600, rH/1000,'color',Colors(2,:),'Linewidth',2,'LineStyle',Style(1));
h(3)=line([0 max(t)/3600], [rGEO/1000 rGEO/1000],...
      'color',Colors(10,:),'Linewidth',2,'LineStyle',Style(2));
h(4)=line([tS(kend)/3600 tS(kend)/3600], [0 rGEO/1000],...
     'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(3));
h(5)=line([tH0 tH0], [0 rGEO/1000],...
     'color',Colors(2,:),'Linewidth',1,'LineStyle',Style(3));
grid on
ylabel('r in 1000 km')
xlim([0 max(t)/3600]); ylim([RE/1000 1.1*rGEO/1000]);
yyaxis right
h(6)= plot(tS/3600, ratioF,'color',Colors(5,:),...
     'Linewidth',1,'LineStyle',Style(1));
hold on
xlabel('t in h')
ylabel('F_S / F_G ')
ylim([ratioF(1) max(ratioF)*1.1]);
legend(h,' Ionenantrieb',' Hohmann-Bahn',' GEO',' t_{ Ionenantrieb}',...
         ' t_{ Hohmann}',' F_S/F_G','location','eastoutside');   
legend box off
ttl=title('Ionenantrieb vs Hohmann-Transfer r(t)');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');


%% Trajektorien-Graphik



figure('Name','Trajektorien')
plot(xS/1000,yS/1000,'color',Colors(3,:),...
     'Linewidth',2,'LineStyle',Style(3));
hold on
plot(xH/1000,yH/1000,'color',Colors(2,:),...
     'Linewidth',2,'LineStyle',Style(3));
plot(r1*xK/1000,r1*yK/1000,'color',Colors(8,:),...
     'Linewidth',1,'LineStyle',Style(1));
plot(rGEO*xK/1000,rGEO*yK/1000,'color',Colors(10,:),...
     'Linewidth',1,'LineStyle',Style(1));
grid on
axis equal
xlim([-rGEO rGEO]*1.2e-3); ylim([-rGEO rGEO]*1.2e-3);
legend(' Ionenantrieb',' Hohmann-Bahn',' LEO', ' GEO');   
legend box off
ttl=title('Trajektorie');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',14);



%% ------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen
 
% DGL
function dY = DGL_IonTriebwerk(t, Y, P1)
% Y(1):x, Y(2):vx, Y(3):y, Y(4):vy
r   = sqrt(Y(1).^2+Y(3).^2);
phi = atan2(Y(3),Y(1));
m   = P1.m0 - P1.mdot.*t;
dY = [Y(2);
      -P1.GME*Y(1)./r^3 - P1.FS*sin(phi)./m ; 
      Y(4);
      -P1.GME*Y(3)./r^3 + P1.FS*cos(phi)./m];  
end 



%Ende Funktionen
% -------------------------------------------------------------------------

