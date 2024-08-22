% -------------------------------------------------------------------------
% BahnParaBestimmung02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Positionen und Bodenspuren 
% eines HEO Satelliten um die Erde aus den Bahnparametern und die
% Beobachtungsdaten einer Bodenstation zu zwei Zeiten. (Vorwärtsrechnung)
% 
% Anschließend wird über das Gauß-Verfahren aus den
% Beobachtungsdaten die Bahnparameter bestimmt. (Rückwärtsrechnung)
%
% Aufgabe mit freundlicher Genehmigung 
% nach Montenbruck, Gill "Satellite Orbits", Springer Verlag 2012
% 
% -------------------------------------------------------------------------

%% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Ausgangs-Parameter

GME = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius

% Beobachtungsstation in Indien
XB = 1344.143;
YB = 6068.601;
ZB = 1429.311;
[lambdaB, phiB, rB]=cart2sph(XB,YB,ZB);
vecXB = [XB;YB;ZB];
lambdaBd = rad2deg(lambdaB);
phiBd    = rad2deg(phiB);


% Einheitsvektoren, Transformationsmatrix Topozentrisch<-->Geozentrisch KOS
eO = [-sin(lambdaB);cos(lambdaB);0];
eN = [-sin(phiB).*cos(lambdaB);-sin(phiB).*sin(lambdaB);cos(phiB)];
eZ = [cos(phiB).*cos(lambdaB);cos(phiB).*sin(lambdaB);sin(phiB)];
U  = [eO,eN,eZ]';
UT = transpose(U);

% Beobachtungszeiten
dt1         = datetime('1999-04-02 00:30:00');
T1          = juliandate(dt1); % Julianisches Datum  ET
dt2         = datetime('1999-04-02 03:00:00');
T2          = juliandate(dt2); % Julianisches Datum  ET
dt3         = datetime('1999-04-02 12:00:00');
T3          = juliandate(dt3); % Julianisches Datum  ET
NPoints1    = 10001;
NPoints3    = 20001;
T_vector1   = linspace(T1,T2,NPoints1);
T_vector2   = linspace(T1,T3,NPoints3);

% Bahnparameter Satellit nach Angabe Montenbruck/Gill

aPn     = 28196.776;
ePn     = 0.7679436;
iPn     = 20.315;
MPn     = 29.236;
omegaPn = 179.425;
OmegaPn = 359.145;
TPn     = 2*pi*sqrt(aPn^3/GME)/86400;

%% Vorwärtsberechnung Berechnung Bahn, Bodenspur und Beobachtungsdaten

% Bodenspur
SatMB.BaPa =[aPn, ePn, iPn, MPn, OmegaPn, omegaPn];
SatMB.Name  = 'HEO Satellit';

SatData  = SatPQR(T_vector1, SatMB, GME);
lat1 (:) = SatData.el;
lat1d(:) = rad2deg(lat1);
theta1d  = GMSTsat(T_vector1);
theta1   = deg2rad(theta1d);
lon1(:)  = SatData.az-theta1;
lon1d(:) = wrapTo360(rad2deg(lon1));

SatData2    = SatPQR(T_vector2, SatMB, GME);
lat2(:)     = SatData2.el;
lat2d(:)    = rad2deg(lat2);
theta2d     = GMSTsat(T_vector2);
theta2      = deg2rad(theta2d);
lon2(:)     = SatData2.az-theta2;
lon2d(:)    = wrapTo360(rad2deg(lon2));



%% Graphische Darstellung
figure('name',SatMB.Name);
gx = geoaxes;
hold on
geoplot(gx, lat2d(:),lon2d(:),'+',...%'MarkerIndices',HourLabel(1).Time(:),...
        'MarkerSize',1,'Color', Colors(9,:),'LineWidth', 1);
geoplot(gx, lat1d(:),lon1d(:),'+',...%'MarkerIndices',HourLabel(1).Time(:),...
        'MarkerSize',1,'Color', Colors(4,:),'LineWidth', 2);
geoplot(gx, phiBd, lambdaBd,'d',...%'MarkerIndices',HourLabel(1).Time(:),...
        'MarkerSize',6,'Color', Colors(3,:),'LineWidth', 2);

geobasemap(gx,'bluegreen')
geolimits([-30 30],[0 +180])
ht = title(strcat(SatMB.Name,' Bodenspur'));
set(ht,'fontweight','normal','fontsize',12)
set(gca,'FontSize',14);

fprintf('\n');
fprintf(SatMB.Name);
fprintf(' \n \n');
fprintf('\n Bahnparameter:');
fprintf('\n a     = %9.3f km  (Große Halbachse)', aPn);
fprintf('\n e     = %9.7f     (Exzentrizität)', ePn);
fprintf('\n i     = %9.3f °   (Inklination)', iPn);
fprintf('\n Omega = %9.3f °   (Rektasz. aufst. Knoten)', OmegaPn);
fprintf('\n omega = %9.3f °   (Argument Perigäum)', omegaPn);
fprintf('\n M     = %9.3f °   (Mittlere Anomalie)', MPn);
fprintf('\n TP    = %9.3f min (Umlaufzeit)', TPn*24*60);
fprintf('\n TP    = %9.3f h   (Umlaufzeit)', TPn*24);
fprintf('\n');
fprintf('\n');


% Geozentrische Satellitenposition zu den Zeiten t1, t2
xA = SatData.xyz(1,1);
yA = SatData.xyz(2,1);
zA = SatData.xyz(3,1);
xE = SatData.xyz(1,end);
yE = SatData.xyz(2,end);
zE = SatData.xyz(3,end);

% Geozentrische Ortsvektoren
vecx1  = [xA;yA;zA];
vecx2  = [xE;yE;zE];
r1     = norm(vecx1);
r2     = norm(vecx2);

% Koordinatentransformation
xS1   = R_z(theta1(1))*vecx1;
xS2   = R_z(theta1(end))*vecx2;
XS1 = xS1 - vecXB; 
XS2 = xS2 - vecXB; 
distance1 = vecnorm(XS1);
distance2 = vecnorm(XS2);
S1  = mtimes(U,XS1);
S2  = mtimes(U,XS2);

% Topozentrische Koordinaten
Az1    = atan2d(S1(1,:),S1(2,:));
hoehe1 = atan2d(S1(3,:),(sqrt(S1(1,:).^2+S1(2,:).^2)));
Az2    = atan2d(S2(1,:),S2(2,:));
hoehe2 = atan2d(S2(3,:),(sqrt(S2(1,:).^2+S2(2,:).^2)));

% Print Ausgabe
fprintf('\n Topozentrische Beobachtung Station Indien: ')
fprintf('\n Länge            = %+6.2f°',  lambdaBd);
fprintf('\n Breite           = %+6.2f°',  phiBd);
fprintf('\n ');
fprintf('\n Beobachtung für %s:\n', dt1);
fprintf('\n Kart. Koordinaten  x1  = [%+10.3f  %+10.3f  %+10.3f] km ', xA, yA, zA);
fprintf('\n Azimut           = %+6.2f°', Az1);
fprintf('\n Höhe             = %+6.2f° ', hoehe1);
fprintf('\n Abstand          = %+6.2f  km', distance1);
fprintf('\n ');
fprintf('\n Beobachtung für %s:\n', dt2);
fprintf('\n Kart. Koordinaten  x2  = [%+10.3f  %+10.3f  %+10.3f] km ', xE, yE, zE);
fprintf('\n Azimut           = %+6.2f°', Az2);
fprintf('\n Höhe             = %+6.2f° ', hoehe2);
fprintf('\n Abstand          = %+6.2f  km', distance2);
fprintf('\n');


%% Rückwärtsrechnung nach Gauss

% Topozentrische Beobachtungdaten Indien
Az1 = deg2rad(132.67);
h1  = deg2rad(32.44);
r1  = 16945.45;
Az2 = deg2rad(123.08);
h2  = deg2rad(50.06);
r2  = 37350.34;

% Berechnung der kartesischen, geozentrischen Koordinaten aus 
% topozentrischen Beobachtungsdaten
[RO, RN, RZ] = sph2cart(pi/2-Az1,h1,r1); %Achtung x-Achse = Ost, 90° geg AZ
vecS         = [RO; RN; RZ];
thetd        = GMSTsat(T1);
thet         = deg2rad(thetd);
vecx1M       = mtimes(transpose(R_z(thet)),(mtimes(UT,vecS)+vecXB));
r1           = norm(vecx1M);
[RO, RN, RZ] = sph2cart(pi/2-Az2,h2,r2);
vecS         = [RO; RN; RZ];
thetd        = GMSTsat(T2);
thet         = deg2rad(thetd);
vecx2M       = mtimes(transpose(R_z(thet)),(mtimes(UT,vecS)+vecXB));
r2           = norm(vecx2M);

vece1 = vecx1M/r1;
vecr0 = vecx2M-dot(vecx2M,vece1)*vece1;
vece0 = vecr0/norm(vecr0);
r0    = norm(vecr0);
W     = cross(vece1,vece0);

inc   = atan2d(sqrt(W(1)^2+W(2)^2),W(3));
Omega = wrapTo360(atan2d(W(1),-W(2)));
u1    = atan2d(vecx1(3),(-vecx1(1)*W(2)+vecx1(2)*W(1)));

tau   = sqrt(GME)*(T2-T1)*86400;
NN    = sqrt(2*(r1*r2+dot(vecx1,vecx2)));
mpara = tau^2/NN^3;
kpara = (r1+r2)/2/NN -0.5;

% ------------------------------------------------------------------------
% Ansatz
% Nullstelle der nichtlinearen Funktion funF(eta) in eta mit 
% der Routine fzero

Dups  = acosd(dot(vecx1,vecx2)/r1/r2);
myfunF = @(eta,P1) 1- eta + P1.m./eta^2*(4*asin(sqrt(P1.m/eta^2-P1.k))- ...
         sin(4*asin(sqrt(P1.m/eta^2-P1.k))))/...
         sin(2*asin(sqrt(P1.m/eta^2-P1.k)))^3;
P1.m = mpara;
P1.k = kpara;
funF = @(x) myfunF(x,P1);  % function of x alone
etaF = fzero(funF,20000);

p     = (etaF*r1*r2*sind(Dups)/tau)^2;
ups1  = atan2d(((p/r1-1)*dot(vecx2M,vece1)/r2-(p/r2-1)),r0*(p/r1-1)/r2);
exz   = (p/r1-1)/cosd(ups1);
a     = p/(1-exz^2);
omega = wrapTo360(u1-ups1);
ups2  = acosd((p/r2-1)/exz);
E1    = atan2(sqrt(1-exz^2)*sind(ups1),(cosd(ups1)+exz));
E2    = atan2(sqrt(1-exz^2)*sind(ups2),(cosd(ups2)+exz));
zetaF = 0.5*(E2-E1);

M     = rad2deg(E1-exz*sin(E1));
TP    = 2*pi*sqrt(a^3/GME)/60;
TPh   = TP/60;

fprintf('\n Bahnparameter bestimmt aus topozentrischer Beobachtung (Indien:')
fprintf(' \n ');
fprintf('\n p     = %9.3f km  (Halbparameter)', p);
fprintf('\n a     = %9.3f km  (Große Halbachse)', a);
fprintf('\n e     = %9.7f     (Exzentrizität)', exz);
fprintf('\n i     = %9.3f °   (Inklination)', inc);
fprintf('\n Omega = %9.3f °   (Rektasz. aufst. Knoten)', Omega);
fprintf('\n omega = %9.3f °   (Argument Perigäum)', omega);
fprintf('\n M     = %9.3f °   (Mittlere Anomalie)', M);
fprintf('\n TP    = %9.3f min (Umlaufzeit)', TP);
fprintf('\n TP    = %9.3f h   (Umlaufzeit)', TPh);
fprintf('\n');
fprintf('\n');

% -------------------------------------------------------------------------
% Alternativer Ansatz I
% Lösung des nichtlinearen Gleichungssystems in eta und zeta mit 
% der Routine fsolve
% Hierzu muß die Vektorfunktion F = func2d(x) der Vektorvariable
% x = [x(1), x(2)] = [eta, zeta] definiert werden.

% ------------------------------------------------------------------------
% Alternativer Ansatz II Iterationslösung für nichtlineare Funktion
 
% eta0  = 12/22+10*sqrt(1+(44/9)*mpara/(kpara+5/6))/22;
% eta1  = eta0+0.1;
% % Erste Iteration
% eta2  = funciter(eta1, eta0, mpara,kpara);
% % Zweite Iteration
% eta0  = eta1;
% eta1  = eta2;
% eta2  = funciter(eta1, eta0, mpara,kpara);
% % Dritte Iteration
% eta0  = eta1;
% eta1  = eta2;
% eta2  = funciter(eta1, eta0, mpara,kpara);

%%  ***************************************************************
%   Funktionen

function F = fun2d(x,P1)
F = zeros(2,1); % Allocate return array
eta  = x(1);
zeta = x(2);
F(1) = P1.m*(2*zeta-sin(2*zeta))/(sin(zeta))^3-eta^2*(eta-1);
F(2) = P1.m/(P1.k+(sin(zeta/2))^2) - eta^2;
end

function eta2 =funciter(eta1, eta0, mpara,kpara)
  f0    = 1-eta0+(mpara/eta0^2)*(4/3+(24/15)*(mpara/eta0^2 - kpara) + ...
        (4*6*8/15/7)*(mpara/eta0^2 - kpara)^2);
  f1    = 1-eta1+(mpara/eta1^2)*(4/3+(24/15)*(mpara/eta1^2 - kpara) + ...
        (4*6*8/15/7)*(mpara/eta1^2 - kpara)^2);
  eta2  = eta1-f1*(eta1-eta0)/(f1-f0);
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
