%%-------------------------------------------------------------------------
% Test.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Bestimmung der Bodenspur der ISS
% und ihrer topozentrischer Beobachtungsdaten von zwei Bodenstationen aus.
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

% Parameter
GME = 398600.4415;  % in km^3/s^2
RE  = 6378;         % in km
dt1 = datetime('2023-03-20 12:00:00');
T1 = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('2023-03-20 13:05:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);


aP   = RE+425.28;
eP   = 0.00185;
iP   = 51.751; 
OmegaP = 47.685;
omegaP = 133;
MPn    =  -51.33;


% Vorwärtsberechnung, Überprüfung
SatHEO.BaPa =[aP , eP, iP , MPn, OmegaP, omegaP];
SatHEO.Name = 'ISS';
SatData = SatPQR(T_vector, SatHEO, GME);
lat1(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon1(:) = wrapTo360(rad2deg(SatData.az)-theta);


aP   = RE+425.28;
eP   = 0.00185;
iP   = 51.751; 
OmegaP = 47.685;
omegaP = 133;
MPn    =  -51.33;


SatHEO.BaPa =[aP , eP, iP , MPn, OmegaP, omegaP];
SatHEO.Name = 'ISS';
SatData = SatPQR(T_vector, SatHEO, GME);
lat2(:) = rad2deg(SatData.el)+5;
theta = GMSTsat(T_vector);
lon2(:) = wrapTo360(rad2deg(SatData.az)-theta);


figure()
gx = geoaxes;
geobasemap(gx,'bluegreen')
geolimits([-50 50],[0 +360])

geoplot(gx, lat1(:),lon1(:),'d',...
       'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);
hold on
geoplot(gx, lat2(:),lon2(:),'o',...
        'MarkerSize',1,'Color', Colors(8,:),'LineWidth', 2);
set(gca,'FontSize',14);
