% -------------------------------------------------------------------------
% Ionantrieb01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Geschwindigkeitsschübe und die Bahnen/Flugzeiten
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
ME  = 5.97e24;       % Masse Erde in kg
ratioTol = 0.025;    % Gültigkeitsbereich für Näherung

% Raumschiff
m0 = 10000;  % Anfangsmasse in kg
FS = 0.6;    % Schub in N
mdot = 3e-3; % Massenausstoß in kg/s

% Hilfsgrößen
phi = linspace(0,2*pi,361);
xK    = cos(phi);
yK    = sin(phi);

%% Berechnung des Geschwindigkeitszuwachs und der Bahn nach Hohmann
v1 = sqrt(GME/r1);
v2 = sqrt(GME/r2);
aH0 = 0.5*(r1+r2);            % Transferellipse Halbachse
vP = sqrt(GME*r2/aH0/r1);     % v bei Perigäum
vA = sqrt(GME*r1/aH0/r2);     % v bei Apogäum
tH0= pi*sqrt(aH0^3/GME)/3600; % Flugzeit für Hohmann Ellipse in h
eccH0 = sqrt(1-r1*r2/aH0^2);  % Exzentrizität

% Geschwindigkeitszuwachs
Dv1 = vP-v1;
Dv2 = v2-vA;
DvHM  = Dv1 + Dv2;

% Bahn der Hohman Ellipse
tH = linspace(0,tH0*3600,1000);
t  = linspace(0,3*tH0*3600,1000);

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

%% Berechnung nach Geschwindigkeitszuwachs kontinuierlicher Schub
v1 = sqrt(GME/r1);
v2 = sqrt(GME/r2);
aH = 0.5*(r1+r2);         % Transferellipse Halbachse
vP = sqrt(GME*r2/aH/r1);  % v bei Perigäum
vA = sqrt(GME*r1/aH/r2);  % v bei Apogäum
% Geschwindigkeitszuwachs
DvIM  = v1 - v2;


%% Berechnung nach Geschwindigkeitszuwachs für Graphik
r1 = RE+h1;
r2 = linspace(RE+h1,RE*100,1000);
v1 = sqrt(GME/r1);
v2 = sqrt(GME./r2);
aH = 0.5*(r1+r2);
vP = sqrt(GME*r2./aH/r1);
vA = sqrt(GME*r1./aH./r2);
Dv1 = (vP/v1-1);
Dv2 = (v2-vA)/v1;
DvH = (Dv1 + Dv2);
DvI = (v1-v2)/v1;
vesc = sqrt(2)-1;
figure('Name','Ionenschub vs Hohmann-Transfer')
semilogx(r2/r1,DvI,'color',Colors(3,:),'Linewidth',2,'LineStyle',Style(1));
hold on
semilogx(r2/r1,DvH,'color',Colors(2,:),'Linewidth',2,'LineStyle',Style(1));
line([rGEO/r1 rGEO/r1], [0 DvIM/v1],'color',Colors(10,:),...
      'linewidth',2, 'Linestyle', ':');
grid on
ylabel('\Deltav / v_1')
xlabel('r_2 / r_1')
legend('\Delta v_I / v_1 (Ionenatrieb)','\Delta v_H / v_1 (Hohmann)',...
       'r_{GEO} = 42157 km  ','location','southeast');   
legend box off
ttl=title('Ionenantrieb vs Hohmann-Transfer - \Delta v');
set(ttl,'FontSize',14,'FontWeight','normal');   
set(gca,'FontSize',14);


%% Graphik Ionenantrieb Abstand Erde und Gravitation/Schub über Zeit
fac1 = sqrt(r1/GME)*FS/mdot;
r     = r1./(1+fac1*log(1-mdot*t/m0)).^2;

kend     = length(t);
kgueltig = length(t);
for k=1:length(t)
   if r(k) >= rGEO
       r(k) = NaN;
       kend = k;
       break
   end
end
for k=kend:length(t)
    r(k) = NaN;
end

rf    = r;
ratioF = FS/m0./(GME./r.^2);

for k=1:kend
  if ratioF(k) >= ratioTol
    r(k)  = NaN;
  else
    rf(k) = NaN;
  end

end
[minr,kgueltig] = max(r);

figure('Name','Ionenantrieb Abstand Erde und Gravitation/Schub über Zeit')
yyaxis left
h(1)=plot(t/3600, r/1000,'color',Colors(3,:),'Linewidth',2,'LineStyle',Style(1));
hold on
plot(t/3600, rf/1000,'color',Colors(3,:),'Linewidth',2,'LineStyle',Style(3));
h(2)=plot(tH/3600, rH/1000,'color',Colors(2,:),'Linewidth',2,'LineStyle',Style(1));
h(3)=line([0 max(t)/3600], [rGEO/1000 rGEO/1000],...
      'color',Colors(10,:),'Linewidth',2,'LineStyle',Style(2));
h(4)=line([t(kend)/3600 t(kend)/3600], [0 rGEO/1000],...
     'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(3));
h(5)=line([tH0 tH0], [0 rGEO/1000],...
     'color',Colors(2,:),'Linewidth',1,'LineStyle',Style(3));
line([0 t(kgueltig)/3600], [r(kgueltig)/1000 r(kgueltig)/1000],...
      'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(3));
grid on
ylabel('r in 1000 km')
xlim([0 max(t)/3600]); ylim([RE/1000 1.1*rGEO/1000]);
% Berechnung Verhältnis Schub/Gravitationskraft
for k=1:kend
   m(k)= m0-mdot*t(k);
end   
for k=kend:length(t)
   m(k)= NaN;
end   
yyaxis right
h(6)= plot(t/3600, ratioF,'color',Colors(5,:),...
     'Linewidth',1,'LineStyle',Style(1));
hold on
line([t(kgueltig)/3600, t(kend)/3600], [ratioTol ratioTol],...
      'color',Colors(5,:),'Linewidth',1,'LineStyle',Style(3));
line([t(kgueltig)/3600 t(kgueltig)/3600], [0.0 ratioTol],...
      'color',Colors(5,:),'Linewidth',1,'LineStyle',Style(3));
xlabel('t in h')
ylabel('F_S / F_G ')
ylim([ratioF(1) ratioTol*2]);
legend(h,' Ionenantrieb',' Hohmann-Bahn',' GEO',' t_{ Ionenantrieb}',...
         ' t_{ Hohmann}',' F_S/F_G','location','eastoutside');   
legend box off
ttl=title('Ionenantrieb vs Hohmann-Transfer r(t)');
set(gca,'FontSize',14);
set(ttl,'FontSize',14,'FontWeight','normal');


%% Trajektorien-Graphik

omegaS= sqrt(GME./r.^3);
xS    = r.*cos(omegaS.*t);
yS    = r.*sin(omegaS.*t);

omegaF= sqrt(GME./rf.^3);
xf    = rf.*cos(omegaF.*t);
yf    = rf.*sin(omegaF.*t);

figure('Name','Trajektorien')
hold on
plot(xS(1:500)/1000,yS(1:500)/1000,'color',Colors(3,:),...
     'Linewidth',2,'LineStyle',Style(3));
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

%% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
