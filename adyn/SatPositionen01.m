% -------------------------------------------------------------------------
% SatPositionen01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Positionen und Bodenspuren 
% von verschiedenen Satelliten um die Erde 
% über die Lösung der Keplergleichung 
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

dt1 = datetime('2000-01-01 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum  ET


%% LEO Satellite
% Berechnung nach Keplerlösung

dt2 = datetime('2000-01-01 09:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);

aPn  = 8000;
iPn  = 55;
ePn  = 0.0;
SatLEO.BaPa =[aPn, ePn, iPn, 0, 0, 90]; %a, e, i, M, Omega, omega
SatLEO.Name = 'LEO';
SatData = SatPQR(T_vector, SatLEO, GME);
lat(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon(:) = wrapTo360(rad2deg(SatData.az)-theta);
figure('name',SatLEO.Name);
gx = geoaxes;
geoplot(gx, lat(:),lon(:),'+',...%'MarkerIndices',HourLabel(1).Time(:),...
        'MarkerSize',1,'Color', Colors(4,:),'LineWidth', 2);
geobasemap(gx,'bluegreen')
geolimits([-60 60],[0 +360])
title(strcat(SatLEO.Name,' Bodenspur'));
set(gca,'FontSize',14);


%% HEO Satellite
%
dt2 = datetime('2000-01-02 16:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);

aPn  = 26555;
SatHEO.BaPa =[aPn, 0.7222, 63.4, 0, 120, 270];
SatHEO.Name = 'HEO';
SatData = SatPQR(T_vector, SatHEO, GME);
lat1(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon1(:) = wrapTo360(rad2deg(SatData.az)-theta);

%% GEO-artiger Satellite (annähernd Geostationärer Satellit)

dt2 = datetime('2000-01-02 12:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);

aPn  = 25000;
SatGEO.BaPa =[aPn, 0.73, 8, 0, 45, 0];
SatGEO.Name = 'GEO';
SatData = SatPQR(T_vector, SatGEO, GME);
lat2(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon2(:) = wrapTo360(rad2deg(SatData.az)-theta);

%% Polar Satellite

dt2 = datetime('2000-01-01 04:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 1001;
T_vector=linspace(T1,T2,NPoints);

aPn  = RE+960;
iPn  = 93;  % in °
OmegaPn = 130;
omegaPn = 0;
ePn     = 0;
SatPol.BaPa =[aPn, ePn, iPn, 0, OmegaPn, omegaPn];
SatPol.Name = 'Polar-Satellit';
SatData = SatPQR(T_vector, SatPol, GME);
lat3(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon3(:) = wrapTo360(rad2deg(SatData.az)-theta);

%% Ausgabe für verschiedene Satelliten-Typen
%

titlestr = 'Bodenspuren verschiedener Satelliten';
figure('name',titlestr);
gx = geoaxes;
geoplot(gx, lat1(:),lon1(:),'+',...
        'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);
text(lat1(1)+2,lon1(1)+2,SatHEO.Name,...
    'VerticalAlignment','bottom','color',Colors(2,:),'FontSize',14);
hold on;
geoplot(gx, lat2(:),lon2(:),'+',...
        'MarkerSize',1,'Color', Colors(3,:),'LineWidth', 2);
text(lat2(1)-25,lon2(1)-20,SatGEO.Name,...
    'VerticalAlignment','bottom','color',Colors(3,:),'FontSize',14);
geoplot(gx, lat3(:),lon3(:),'+',...
        'MarkerSize',1,'Color', Colors(4,:),'LineWidth', 2);
text(lat3(1),lon3(1)+2,SatPol.Name,...
    'VerticalAlignment','bottom','color',Colors(4,:),'FontSize',14);

geolimits([-80 80],[0 +360])
geobasemap(gx,'bluegreen')
title(titlestr);
set(gca,'FontSize',14);

%% Verschiedene GEO-synchrone Satelliten

dt1 = datetime('2000-01-01 00:00:00');
dt2 = datetime('2000-01-02 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum  ET
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 501;
T_vector=linspace(T1,T2,NPoints);
ePn = 0.2;
TStar = 86164; %Sterntag in sec
aPn   = (TStar^2*GME/4/pi^2)^(1/3);
SatGEO.BaPa =[aPn, ePn, -20, 30, 30, 360];
SatGEO.Name1 = 'GEO1';
SatData = SatPQR(T_vector, SatGEO, GME);
latG1(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lonG1(:) = wrapTo180(rad2deg(SatData.az)-theta);

SatGEO.BaPa =[aPn, ePn, -30, 30, -35, 35];
SatGEO.Name2 = 'GEO2';
SatData = SatPQR(T_vector, SatGEO, GME);
latG2(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lonG2(:) = wrapTo180(rad2deg(SatData.az)-theta);

SatGEO.BaPa =[aPn, ePn, 45, 30, -100, 30];
SatGEO.Name3 = 'GEO3';
SatData = SatPQR(T_vector, SatGEO, GME);
latG3(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lonG3(:) = wrapTo180(rad2deg(SatData.az)-theta);

titlestr = 'Bodenspuren verschiedener geosynchroner-Satelliten';
figure('name',titlestr);
gx = geoaxes;
geoplot(gx, latG1(:),lonG1(:),'+',...
        'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);
text(min(latG1),max(lonG1),SatGEO.Name1,...
    'VerticalAlignment','bottom','color',Colors(2,:),'FontSize',14);
hold on;
geoplot(gx, latG2(:),lonG2(:),'+',...
        'MarkerSize',1,'Color', Colors(3,:),'LineWidth', 2);
text(max(latG2),min(lonG2),SatGEO.Name2,...
    'VerticalAlignment','bottom','color',Colors(3,:),'FontSize',14);
geoplot(gx, latG3(:),lonG3(:),'+',...
        'MarkerSize',1,'Color', Colors(4,:),'LineWidth', 2);
text(max(latG3),max(lonG3),SatGEO.Name3,...
    'VerticalAlignment','bottom','color',Colors(4,:),'FontSize',14);

geolimits([-45 45],[-120 +120])
geobasemap(gx,'bluegreen')
title(titlestr);
set(gca,'FontSize',14);




% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
