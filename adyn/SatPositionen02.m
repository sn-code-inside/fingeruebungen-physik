% -------------------------------------------------------------------------
% SatPositionen02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Bahn, die Bodenspur, 
% und die Bahnparameter
% der ISS auf Basis von Daten der NASA
% (https://spotthestation.nasa.gov/trajectory_data.cfm)
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


%% Parameter und Datenimport
GME = 398600.4415;   % G*ME in 
RE  = 6378;          % Erdradius

ISS = ImportISSData('ISS.txt');  % Einlesen der NASA Daten
% ISS state vectors in the Mean of J2000 (J2K) reference frame are 
% listed at four-minute intervals spanning a total length of 15 days. 
% During reboosts (translation maneuvers), the state vectors are reported 
% in two-second intervals. Each state vector lists the time in UTC; 
% position X, Y, and Z in km; and velocity X, Y, and Z in km/s.

dt1 = datetime(ISS.Date(1));        % Startdatum der NASA-Daten
T1 = juliandate(dt1);               % Julianisches Datum  ET

% Koordinaten- und Geschwindigkeitsvektoren aus NASA Daten
x = ISS.x;
y = ISS.y;
z = ISS.z;
vx = ISS.vx;
vy = ISS.vy;
vz = ISS.vz;

%% Berechnung von Az, Alt, Radius aus kartesischen NASA-Koordinaten

[az,el,r] = cart2sph(x,y,z);
hS = r - RE;  % Höhe über mittleren Erdradius
% Hier könnte man auch die Abplattung berücksichtigen. Dann wird hS auch 
% noch eine Funktion von del und nicht nur von r
v  = sqrt(vx.^2+vy.^2+vz.^2);
az    = rad2deg(az);
del   = rad2deg(el);

% Graphik der Variation der Höhe der ISS über die Zeit 
% leichte Exzentrizität der Keplerbahn und Abplattungseffekt der Erde

titlestr = strjoin(["Höhe u. Deklination ISS"; "(NASA-Daten)           "]);
figure('name',titlestr)
yyaxis left
hS  = r - RE;
h(1) = plot(ISS.Date,hS,'LineWidth',2,'LineStyle', Style(3));
% Korrektur Erdabplattung
f   = 1 -0.996647;  % Abplattungskonstante ca 1/300
xE  = RE.*cosd(del).*cosd(az);
yE  = RE.*cosd(del).*sind(az);
zE  = RE.*sind(del);
hSk = r - sqrt(xE.^2+yE.^2+zE.^2/(1-f)^2);
hold on;
ylim([350,450]);
h(2) = plot(ISS.Date,hSk,'LineWidth',2,'LineStyle', Style(1));
ylabel('Höhe über Erde in km');
yyaxis right
h(3) =plot(ISS.Date,del,'LineWidth',2);
line([dt1,max(ISS.Date)],[0,0],'LineWidth',1);
ylabel('Deklination °');
grid on;
ylim([-75,75]);
legend(h,'Höhe über Erdkugel', 'Höhe über Ellipsoid', 'Deklination',...
       'location','south');
legend box off
hp1 = title(titlestr);
set(hp1,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',14);

% %%  3D-Graphik der Bahn
% titlestr = strjoin(["ISS-Bahn (3D) "; "(NASA-Daten)  "]);
% figure('name',titlestr)
% plot3(x,y,z);
% axis equal
% hold on
% set(gca,'FontSize',14);
% maxR = 1.5*RE;
% xlim([-maxR maxR]);
% ylim([-maxR maxR]);
% zlim([-maxR maxR]);
% grid on
% xlabel('x in km');
% ylabel('y in km');
% zlabel('z in km');
% xl = xlim;
% yl = ylim;
% [XB,YB] = meshgrid(xl,yl);
% surf(XB,YB,zeros(size(XB)))
% shading flat
% alpha 0.2
% [xs,ys,zs] = sphere(60);  
% surf(RE*xs,RE*ys,RE*zs) % Erdkugel
% colormap winter;
% plot3(r.*cosd(del).*cosd(az),r.*cosd(del).*sind(az),r.*sind(del));
% nend = 30;
% for k=1:nend  %Einzelne Datenpunkte aus der NASA-Datei
%     plot3(x(k),y(k),z(k),'s','Color','r', 'MarkerEdgeColor','r',...
%           'MarkerSize',5,'Linewidth',4);
%     pause(0.1);
% end
% hp1 = title(titlestr,'FontSize',12);
% set(hp1,'FontSize',14,'FontWeight','normal'); 
% set(gca,'FontSize',14);
% 
%% ISS-Bodenspur aus NASA-Daten

dt2 = datetime(ISS.Date);
T2 = juliandate(dt2); % Julianisches Datum  ET
lat_NASA = del;
theta = GMSTsat(T2);
lon_NASA = wrapTo360(az-theta);

% Graphische Ausgabe Bodenspuren NASA
titlestr = strjoin(["ISS-Bodenspur "; "(NASA-Daten und Kepler-Berechnung)  "]);
figure('name',titlestr);
gx = geoaxes;
nend  = 30;
equat = linspace(0,360,37);
geoplot(gx, lat_NASA(1:nend),lon_NASA(1:nend),...
        '+','Color','r', 'MarkerEdgeColor','r',...
          'MarkerSize',6,'Linewidth',2);
text(lat_NASA(1:nend),lon_NASA(1:nend)+1, string(dt2(1:nend),...
    'HH:mm'),'FontSize',10,'Color',Colors(4,:));
for k=2:nend
  if lon_NASA(k) > lon_NASA(k-1)
  line([lat_NASA(k) lat_NASA(k-1)],[lon_NASA(k) lon_NASA(k-1)],...
       'Color','r');
  end
end
hold on
geoplot(gx, 0*equat, equat, 'color','b','LineWidth',1);
geobasemap(gx,'bluegreen')
geolimits([-50 50],[0 +360])
text(-65,190, strcat(string(dt1,'yyyy-MM-dd'),' (UTC)'),...
     'FontSize',12,'Color',Colors(4,:));
hp1 = title(titlestr,'FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',14);


%% Bahnparameterberechnung aus Messungen 

NSample = 1;           %Nr des Messwertes aus dem Datensatz
x = ISS.x(NSample);
y = ISS.y(NSample);
z = ISS.z(NSample);
vx = ISS.vx(NSample);
vy = ISS.vy(NSample);
vz = ISS.vz(NSample);
r  = sqrt(x^2+y^2+z^2);
hS = r - RE;
v  = sqrt(vx^2+vy^2+vz^2);

%Drehimpuls
Lx =y*vz-z*vy;
Ly =z*vx-x*vz;
Lz =x*vy-y*vx;

L = sqrt(Lx^2+Ly^2+Lz^2);

Wx = Lx/L;
Wy = Ly/L;
Wz = Lz/L;

iP     = atan(sqrt(Wx^2+Wy^2)/Wz);
OmegaP = atan2(Wx,-Wy);

p      = L^2/GME;
aP     = 1/(2/r - v^2/GME);
TPn    = 2*pi*sqrt(aP^3/GME);
TPmin  = 2*pi*sqrt(aP^3/GME)/60;
nP     = sqrt(GME/aP^3);
eP     = sqrt(1-p/aP);
E      = atan2((x*vx+y*vy+z*vz)/aP^2/nP,(1-r/aP));
M      = E - eP*sin(E);
uP     = atan2(z,(-x*Wy+y*Wx));
wA     = atan2(sqrt(1-eP^2)*sin(E),(cos(E)-eP));
omegaP = uP - wA;
iP     = rad2deg(iP);
OmegaP = wrapTo360(rad2deg(OmegaP));
omegaP = wrapTo360(rad2deg(omegaP));
MP     = wrapTo360(rad2deg(M));
% Vorw�rtsberechnung, �berpr�fung
% Parameter
GME = 398600.4415;  % in km^3/s^2
RE  = 6378;         % in km
dt1 = datetime('2023-03-20 12:00:00');
T1 = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('2023-03-20 13:56:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);
% Parameter
aP     = RE+425.28;
eP     = 0.00185;
iP     = 51.751; 
OmegaP = 47.685;
omegaP = 124;
%!!! Mittlere Anomalie wird angepasst, damit Anfangswert mit NASA Daten
%�bereinstimmt. Siehe Bemerkung am Ende der �bung!
MP     = wrapTo360(315.75);  
% Vorw�rtsberechnung
SatHEO.BaPa =[aP , eP, iP , MP, OmegaP, omegaP];
SatHEO.Name = 'ISS';
SatData = SatPQR(T_vector, SatHEO, GME);
lat1(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon1(:) = wrapTo360(rad2deg(SatData.az)-theta);
geoplot(gx, lat1(:),lon1(:),'d',...
        'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);

    
% Asudruck
fprintf('\n');
fprintf('Bahnparameter ISS\n');
fprintf('\n t = %s', dt1);
fprintf('\n r = %8.2f km', r);
fprintf('\n v = %8.2f km/s', v);
fprintf('\n h = %8.2f km', hS);
fprintf('\n');
fprintf('\n a     = %8.2f km  (Große Halbachse)', aP);
fprintf('\n e     = %8.5f     (Exzentrizität)', eP);
fprintf('\n i     = %8.2f     (Inklination)', iP);
fprintf('\n Omega = %8.2f     (Rektasz. aufst. Knoten)', OmegaP);
fprintf('\n omega = %8.2f     (Argument Perigäum)', omegaP);
fprintf('\n M     = %8.2f     (Mittlere Anomalie)', MP);
fprintf('\n TP    = %8.2f min (Umlaufzeit)', TPmin);
fprintf('\n');
fprintf('\n');
fprintf('\n Achtung !!! Auf Grund der nahezu kreisförmigen Bahn');
fprintf('\n k�nnen einige Bahnparameter, insbesondere das Argument');
fprintf('\n des Perigäums und die Exzentrizität, nicht genau definiert');
fprintf('\n bzw. nur sehr ungenau zu bestimmen sein. ');
fprintf('\n Man muss dann auf Spezialverfahren zurück greifen. ');
fprintf('\n');
fprintf('\n Siehe Kap. 3.3.4.2 Bestimmung Bahnparameter bei kreisf�rmigen');
fprintf('\n und �quatorialen BahnenKapitel und Zusatzaufgabe.) ');
fprintf('\n');

   

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

