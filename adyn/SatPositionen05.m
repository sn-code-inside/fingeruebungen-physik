%%-------------------------------------------------------------------------
% SatPositionen03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Bestimmung der topozentrische Beobachtungsdaten eines polaren
% Satelliten von zwei Bodenstationen aus.
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
dt1 = datetime('2020-03-11 16:00:0');
T1 = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('2020-03-12 04:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);


% Stationen

phiB1    = 10;  %1 Indien  
lambdaB1 = 75;  
phiB2    = 55;  %2 Russland
lambdaB2 = 75;  
LocStr(2,:) = "Russland";
LocStr(1,:) = "Indien";

%% Vorwärtsrechnung aus Bahnparametern

aP   = 26555;
eP   = 0.72;
iP   = 63.4; 
OmegaP = 135;
omegaP = 270;
TPn  = 2*pi*sqrt(aP^3/GME)/86400;
M    = 0;

SatHEO.BaPa =[aP , eP, iP , M, OmegaP, omegaP];
SatHEO.Name = 'Molniya-Satellit';
SatData = SatPQR(T_vector, SatHEO, GME);
lat(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon(:) = wrapTo360(rad2deg(SatData.az)-theta);

fprintf('\n');
fprintf('Bahnparameter für Bahnberechnung nach Kepler\n');
fprintf('\n Zeit  = %s - %s ', dt1, dt2);
fprintf('\n');
fprintf('\n a     = %8.2f km  (Große Halbachse)', aP);
fprintf('\n e     = %8.5f     (Exzentrizität)', eP);
fprintf('\n i     = %8.2f °   (Inklination)', iP);
fprintf('\n Omega = %8.2f °   (Rektasz. aufst. Knoten)', OmegaP);
fprintf('\n omega = %8.2f °   (Argument Perigäum)', omegaP);
fprintf('\n M     = %8.2f °   (Mittlere Anomalie)', M);
fprintf('\n TP    = %8.2f min (Umlaufzeit)', TPn*24*60);
fprintf('\n');
fprintf('\n');

%% Graphik Ground Track

titlestr = strcat(SatHEO.Name,' Ground Track');
figure('name',titlestr);
gx = geoaxes;
geoplot(gx, lat(:),lon(:),'+',...
        'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);
hold on
geoplot(gx, phiB1,lambdaB1,'o',...
        'MarkerSize',5,'Color', Colors(4,:),'LineWidth', 2);
geoplot(gx, phiB2,lambdaB2,'d',...
        'MarkerSize',5,'Color', Colors(4,:),'LineWidth', 2);
geobasemap(gx,'bluegreen')
geolimits([-60 60],[-120 +240])
title(titlestr);
set(gca,'FontSize',14);

%% Topozentrische Beobachtungsberechnung Indien

r1     = SatData.r;
xS1 = r1.*[cosd(lat).*cosd(lon);cosd(lat).*sind(lon);sind(lat)];


xB1 = RE*[cosd(phiB1).*cosd(lambdaB1);cosd(phiB1).*sind(lambdaB1);sind(phiB1)];
XS1 = xS1- xB1; 
distance1 = vecnorm(XS1);

eO = [-sind(lambdaB1);cosd(lambdaB1);0];
eN = [-sind(phiB1).*cosd(lambdaB1);-sind(phiB1).*sind(lambdaB1);cosd(phiB1)];
eZ = [cosd(phiB1).*cosd(lambdaB1);cosd(phiB1).*sind(lambdaB1);sind(phiB1)];
U  = [eO,eN,eZ]';
S  = mtimes(U,XS1);
Az1    = atan2d(S(1,:),S(2,:));
hoehe1 = atan2d(S(3,:),(sqrt(S(1,:).^2+S(2,:).^2)));


%% Topozentrische Beobachtungsberechnung Russland

xB2 = RE*[cosd(phiB2).*cosd(lambdaB2);cosd(phiB2).*sind(lambdaB2);sind(phiB2)];
XS2 = xS1- xB2; 
distance2 = vecnorm(XS2);

eO = [-sind(lambdaB2);cosd(lambdaB2);0];
eN = [-sind(phiB2).*cosd(lambdaB2);-sind(phiB2).*sind(lambdaB2);cosd(phiB2)];
eZ = [cosd(phiB2).*cosd(lambdaB2);cosd(phiB2).*sind(lambdaB2);sind(phiB2)];
U  = [eO,eN,eZ]';
S  = mtimes(U,XS2);
Az2    = atan2d(S(1,:),S(2,:));
hoehe2 = atan2d(S(3,:),(sqrt(S(1,:).^2+S(2,:).^2)));


for k = 1:length(Az1)
   if hoehe1(k)< 0 
       Az1(k) = NaN;
       distance1(k) = NaN;
   end
   if hoehe2(k)< 0 
       Az2(k) = NaN;
       distance2(k) = NaN;
   end
end

%% Graphische Ausgabe

titlestr = strcat(SatHEO.Name,' Topozentr. Beobachtung');
figure('name',titlestr);
subplot(2,1,1)
yyaxis left
plot(datetime(SatData.Time,'ConvertFrom','juliandate'),hoehe1,...
    'LineWidth',2);
hold on
plot(datetime(SatData.Time,'ConvertFrom','juliandate'),hoehe2,...
    'LineWidth',2,'LineStyle',Style(3));
ylabel('Hoehe in °')
ylim([0,90]);
yyaxis right
ylabel('Azimut in °')
p(1)=plot(datetime(SatData.Time,'ConvertFrom','juliandate'),Az1,...
    'LineWidth',2);
hold on
p(2)=plot(datetime(SatData.Time,'ConvertFrom','juliandate'),Az2,...
    'LineWidth',2,'LineStyle',Style(3));
grid on;
hp1 = title('Höhe und Azimut','FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:2),LocStr(1:2,:),'location','south','NumColumns',1);
set(hp2,'FontSize',14,'FontWeight','normal');
legend box off
set(gca,'FontSize',14);


subplot(2,1,2)
plot(datetime(SatData.Time,'ConvertFrom','juliandate'),...
     distance1,'color',Colors(3,:),'LineWidth',2);
hold on
plot(datetime(SatData.Time,'ConvertFrom','juliandate'),...
     distance2,'color',Colors(3,:),'LineWidth',2,'LineStyle',Style(3));
ylabel('Entfernung in km')
ylim([0,1.1*max(max(distance1),max(distance2))]);
grid on;
hp1 = title('Entfernung','FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(LocStr(1:2,:),'location','south','NumColumns',1);
set(hp2,'FontSize',14,'FontWeight','normal');
legend box off
set(gca,'FontSize',14);


%% Ende Programm


