% -------------------------------------------------------------------------
% SatPositionenTestMontenbruck.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Positionen und Bodenspuren 
% von verschiedenen Satelliten um die Erde 
% Ã¼ber die LÃ¶sung der Keplergleichung 
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
RE  = 6378;          % Erdradius




%% Montenbruck Satellite
%
dt1 = datetime('1999-04-02 00:30:00');
T1 = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('1999-04-02 03:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 1001;
T_vector=linspace(T1,T2,NPoints);
Timevector = datetime(T_vector,'ConvertFrom','juliandate');
aPn  = 28196.776;
ePn  = 0.7679436;
iPn  = 20.315;
MPn  = 29.236;
omegaPn = 179.425;
OmegaPn = 359.145;
TPn  = 2*pi*sqrt(aPn^3/GME)/86400;

SatMB.BaPa =[aPn, ePn, iPn, MPn, OmegaPn, omegaPn];
SatMB.Name = 'Montenbruck';
SatData = SatPQR_mod(T_vector, SatMB, GME);
lat(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon(:) = wrapTo360(rad2deg(SatData.az)-theta);
figure('name',SatMB.Name);
gx = geoaxes;
geoplot(gx, lat(:),lon(:),'+',...%'MarkerIndices',HourLabel(1).Time(:),...
        'MarkerSize',1,'Color', Colors(4,:),'LineWidth', 2);
geobasemap(gx,'bluegreen')
geolimits([-60 60],[0 +360])
title(strcat(SatMB.Name,' Bodenspur'));
set(gca,'FontSize',14);

xP = SatData.xyz(1,1);
yP = SatData.xyz(2,1);
zP = SatData.xyz(3,1);

xE = SatData.xyz(1,end);
yE = SatData.xyz(2,end);
zE = SatData.xyz(3,end);

fprintf('\n');
fprintf(strcat(SatMB.Name,'-Satellit: Bahnparameter'),' \n \n');
fprintf('\n \n Zeit  = %s', dt1);
fprintf('\n');
fprintf('\n a     = %9.3f km  (Große Halbachse)', aPn);
fprintf('\n e     = %9.7f     (Exzentrizität)', ePn);
fprintf('\n i     = %9.3f °   (Inklination)', iPn);
fprintf('\n Omega = %9.3f °   (Rektasz. aufst. Knoten)', OmegaPn);
fprintf('\n omega = %9.3f °   (Argument Perigäum)', omegaPn);
fprintf('\n M     = %9.3f °   (Mittlere Anomalie)', MPn);
fprintf('\n TP    = %9.3f min (Umlaufzeit)', TPn*24*60);
fprintf('\n TP    = %9.3f h   (Umlaufzeit)', TPn*24);
fprintf('\n');
fprintf('\n Kart. Koordinaten  xv0  = [%+10.3f  %+10.3f  %+10.3f] km ', xP, yP, zP);
fprintf('\n Kart. Koordinaten  xvE  = [%+10.3f  %+10.3f  %+10.3f] km ', xE, yE, zE);
fprintf('\n Abstand Erdzentrum r    =  %+10.3f  km', norm([xP yP zP]));
fprintf('\n Abstand Erdzentrum r    =  %+10.3f  km', SatData.r(1));
fprintf('\n');
fprintf('\n');

% Beobachtungsstation Indien
XB = 1344.143;
YB = 6068.601;
ZB = 1429.311;
[lambdaB, phiB, rB]=cart2sph(XB,YB,ZB);
phiB    = rad2deg(phiB);
lambdaB = rad2deg(lambdaB);



%% Topozentrische Beobachtungsberechnung Indien

vecx1  = [xP;yP;zP];
vecx1  = [xE;yE;zE];
r1     = norm(vecx1);
xS1    = r1*[cosd(lat(1)).*cosd(lon(1));cosd(lat(1)).*sind(lon(1));sind(lat(1))];
xS1    = r1*[cosd(lat(end)).*cosd(lon(end));cosd(lat(end)).*sind(lon(end));sind(lat(end))];

xB  = [XB;YB;ZB];
XS1 = xS1 - xB; 
distance1 = vecnorm(XS1);

eO = [-sind(lambdaB);cosd(lambdaB);0];
eN = [-sind(phiB).*cosd(lambdaB);-sind(phiB).*sind(lambdaB);cosd(phiB)];
eZ = [cosd(phiB).*cosd(lambdaB);cosd(phiB).*sind(lambdaB);sind(phiB)];
U  = [eO,eN,eZ]';
S  = mtimes(U,XS1);
Az1    = atan2d(S(1,:),S(2,:));
hoehe1 = atan2d(S(3,:),(sqrt(S(1,:).^2+S(2,:).^2)));

fprintf('\n Station Indien: ')
fprintf('\n Länge            = %+6.2f°',  lambdaB);
fprintf('\n Breite           = %+6.2f°',  phiB);
fprintf('\n Azimut           = %+6.2f°', Az1);
fprintf('\n Höhe             = %+6.2f° ', hoehe1);
fprintf('\n Abstand          = %+6.2f  km', distance1);
fprintf('\n');
fprintf('\n');


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
