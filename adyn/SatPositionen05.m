%%-------------------------------------------------------------------------
% SatPositionen03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Finger√ºbungen" von Michael Kaschke und Holger Cartarius
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
dt1 = datetime('2020-03-11 00:00:0');
T1 = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('2020-03-11 16:00:00');
T2 = juliandate(dt2); % Julianisches Datum  ET
NPoints = 5001;
T_vector=linspace(T1,T2,NPoints);

xB = 1344;
yB = 6068;
zB = 1429;

% Stationen
phi2     = atan2d(zB,sqrt(xB^2+yB^2));
lambda2  = atan2d(yB,xB);

phiB     = [51,phi2];  %1 Gˆrlitz  2 %Indien
lambdaB  = [15,lambda2];  
LocStr(1,:) = "Gˆrlitz (D)";
LocStr(2,:) = "Indien";
%% Vorw√§rtsrechnung aus Bahnparametern

% Berechnung Bahn nach Keplerl√∂sung
aP   = 28000;
eP   = 0.75;
iP   = 21; 
OmegaP = 360;
omegaP = 180;
TPn  = 2*pi*sqrt(aP^3/GME)/86400;
M    = 20;
SatHEO.BaPa =[aP , eP, iP , M, OmegaP, omegaP];
SatHEO.Name = 'HEO-Satellit';
SatData = SatPQR(T_vector, SatHEO, GME);
lat1(:) = rad2deg(SatData.el);
theta = GMSTsat(T_vector);
lon1(:) = wrapTo360(rad2deg(SatData.az)-theta);

fprintf('\n');
fprintf('Bahnparameter f√ºr Bahnberechnung nach Kepler\n');
fprintf('\n Zeit  = %s', dt1);
fprintf('\n');
fprintf('\n a     = %8.2f km  (Gro√üe Halbachse)', aP);
fprintf('\n e     = %8.5f     (Exzentrizit√§t)', eP);
fprintf('\n i     = %8.2f ∞   (Inklination)', iP);
fprintf('\n Omega = %8.2f ∞    (Rektasz. aufst. Knoten)', OmegaP);
fprintf('\n omega = %8.2f ∞   (Argument Perig√§um)', omegaP);
fprintf('\n M     = %8.2f ∞   (Mittlere Anomalie)', M);
fprintf('\n TP    = %8.2f min (Umlaufzeit)', TPn*24*60);
fprintf('\n');
fprintf('\n');

%% Graphik Ground Track

titlestr = strcat(SatHEO.Name,' Ground Track');
figure('name',titlestr);
gx = geoaxes;
geoplot(gx, lat1(:),lon1(:),'+',...
        'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);
hold on
geoplot(gx, phiB(1),lambdaB(1),'o',...
        'MarkerSize',5,'Color', Colors(4,:),'LineWidth', 2);
geoplot(gx, phiB(2),lambdaB(2),'d',...
        'MarkerSize',5,'Color', Colors(4,:),'LineWidth', 2);
geobasemap(gx,'bluegreen')
geolimits([-60 60],[0 +360])
title(titlestr);
set(gca,'FontSize',14);

%% Topozentrische Beobachtungsberechnung

xS  = aP*[cosd(lat1).*cosd(lon1);cosd(lat1).*sind(lon1);sind(lat1)];
xB = RE*([cosd(phiB).*cosd(lambdaB);cosd(phiB).*sind(lambdaB);sind(phiB)]);

%Gˆrlitz
XS1 = xS - xB(:,1); 
distance1 = vecnorm(XS1);
eO = [-sind(lambdaB(1));cosd(lambdaB(1));0];
eN = [-sind(phiB(1)).*cosd(lambdaB(1));-sind(phiB(1)).*sind(lambdaB(1));cosd(phiB(1))];
eZ = [cosd(phiB(1)).*cosd(lambdaB(1));cosd(phiB(1)).*sind(lambdaB(1));sind(phiB(1))];
U  = [eO,eN,eZ]';
S  = mtimes(U,XS1);
Az1    = atan2d(S(1,:),S(2,:));
hoehe1 = atan2d(S(3,:),(sqrt(S(1,:).^2+S(2,:).^2)));

%Umea
XS2 = xS - xB(:,2); 
distance2 = vecnorm(XS2);
eO = [-sind(lambdaB(2));cosd(lambdaB(2));0];
eN = [-sind(phiB(2))*cosd(lambdaB(2));-sind(phiB(2))*sind(lambdaB(2));cosd(phiB(2))];
eZ = [cosd(phiB(2))*cosd(lambdaB(2));cosd(phiB(2))*sind(lambdaB(2));sind(phiB(2))];
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
ylabel('Hoehe in ∞')
ylim([0,90]);
yyaxis right
ylabel('Azimut in ∞')
p(1)=plot(datetime(SatData.Time,'ConvertFrom','juliandate'),Az1,...
    'LineWidth',2);
hold on
p(2)=plot(datetime(SatData.Time,'ConvertFrom','juliandate'),Az2,...
    'LineWidth',2,'LineStyle',Style(3));
grid on;
hp1 = title('Hˆhe und Azimut','FontSize',12);
set(hp1,'FontSize',14,'FontWeight','normal'); 
hp2=legend(p(1:2),LocStr(1:2,:),'location','northeast','NumColumns',1);
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
hp2=legend(LocStr(1:2,:),'location','southeast','NumColumns',1);
set(hp2,'FontSize',14,'FontWeight','normal');
legend box off
set(gca,'FontSize',14);


%% Ende Programm


