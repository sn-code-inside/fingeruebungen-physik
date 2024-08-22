%%-------------------------------------------------------------------------
% SatPositionen03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Bestimmung der topozentrische Beobachtungsdaten eines MEO
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

% Stationen
Breite      = [51,60];  %1 Görlitz  2 %Marienhamn (Schweden)
Laenge      = [15,20];  
LocStr(1,:) = "Görlitz (D)";
LocStr(2,:) = "Marienham (S)";

GME = 398600.4415;   % G*ME in km^3/s^2
RE  = 6378;          % Erdradius


% Stationen
Breite      = [51,60];  %1 Görlitz  2 %Marienhamn (Schweden)
Laenge      = [15,20];  
LocStr(1,:) = "Görlitz D";
LocStr(2,:) = "Marienham S";

% Beobachtungsstation Görlitz
phiBd        = Breite(1);
phiB         = deg2rad(phiBd);
lambdaBd     = Laenge(1);
lambdaB      = deg2rad(lambdaBd);
[XB, YB, ZB] = sph2cart(lambdaB,phiB,RE);
vecXB        = [XB;YB;ZB];

% % Beobachtungsstation Marienham
% phiBd        = Breite(2);
% phiB         = deg2rad(phiBd);
% lambdaBd     = Laenge(2);
% lambdaB      = deg2rad(lambdaBd);
% [XB, YB, ZB] = sph2cart(lambdaB,phiB,RE);
% vecXB        = [XB;YB;ZB];


% Einheitsvektoren, Transformationsmatrix Topozentrisch<-->Geozentrisch KOS
eO = [-sin(lambdaB);cos(lambdaB);0];
eN = [-sin(phiB).*cos(lambdaB);-sin(phiB).*sin(lambdaB);cos(phiB)];
eZ = [cos(phiB).*cos(lambdaB);cos(phiB).*sin(lambdaB);sin(phiB)];
U  = [eO,eN,eZ]';
UT = transpose(U);


% Beobachtungszeiten
dt1 = datetime(['2020-01-01 03:16:00']);
T1  = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('2020-01-01 03:25:00');
T2  = juliandate(dt2); % Julianisches Datum  ET
dt3 = datetime('2020-01-01 04:30:00');
T3  = juliandate(dt3); % Julianisches Datum  ET
NPoints1     = 1001;
NPoints2    = 10001;
T_vector1   = linspace(T1,T2,NPoints1);
T_vector2   = linspace(T1,T3,NPoints2);

%% Bahnparameter Satellit(Vorwärtsrechnung)

aP     = RE+1500;
eP     = 0.05;
iP     = 65; 
OmegaP = 135;
omegaP = 115;
TPn    = 2*pi*sqrt(aP^3/GME)/86400;
M      = 290;
SatHEO.BaPa = [aP , eP, iP , M, OmegaP, omegaP];
SatHEO.Name = 'MEO-Satellit';

SatData1    = SatPQR(T_vector1, SatHEO, GME);
lat1(:)     = SatData1.el;
lat1d(:)    = rad2deg(lat1);
theta1d     = GMSTsat(T_vector1);
theta1      = deg2rad(theta1d);
lon1(:)     = SatData1.az-theta1;
lon1d(:)    = wrapTo360(rad2deg(lon1));

SatData2    = SatPQR(T_vector2, SatHEO, GME);
lat2(:)     = SatData2.el;
lat2d(:)    = rad2deg(lat2);
theta2d     = GMSTsat(T_vector2);
theta2      = deg2rad(theta2d);
lon2(:)     = SatData2.az-theta2;
lon2d(:)    = wrapTo360(rad2deg(lon2));


fprintf('\n');
fprintf('Bahnparameter fÃ¼r Bahnberechnung nach Kepler\n');
fprintf('\n Zeit  = %s', dt1);
fprintf('\n');
fprintf('\n a     = %8.2f km  (GroÃŸe Halbachse)', aP);
fprintf('\n e     = %8.5f     (ExzentrizitÃ¤t)', eP);
fprintf('\n i     = %8.2f     (Inklination)', iP);
fprintf('\n Omega = %8.2f     (Rektasz. aufst. Knoten)', OmegaP);
fprintf('\n omega = %8.2f     (Argument PerigÃ¤um)', omegaP);
fprintf('\n TP    = %8.2f min (Umlaufzeit)', TPn*24*60);
fprintf('\n');
fprintf('\n');

%% Graphik Ground Track

titlestr = strcat(SatHEO.Name,' Ground Track');
figure('name',titlestr);
gx = geoaxes;
hold on
geoplot(gx, lat2d(:),lon2d(:),'+',...
        'MarkerSize',1,'Color', Colors(4,:),'LineWidth', 2);
geoplot(gx, lat1d(:),lon1d(:),'+',...
        'MarkerSize',1,'Color', Colors(2,:),'LineWidth', 2);
geoplot(gx, rad2deg(phiB),rad2deg(lambdaB),'o',...
        'MarkerSize',5,'Color', Colors(4,:),'LineWidth', 2);
geobasemap(gx,'bluegreen')
geolimits([-90 90],[0 +360])
title(titlestr);
set(gca,'FontSize',14);


%% Berechnung topozentrischer Beobachtungsdaten

% Geozentrische Satellitenposition zu den Zeiten t1, t2
xA = SatData1.xyz(1,1);
yA = SatData1.xyz(2,1);
zA = SatData1.xyz(3,1);
xE = SatData1.xyz(1,end);
yE = SatData1.xyz(2,end);
zE = SatData1.xyz(3,end);

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



%% Graphische Ausgabe

xS  = aP*[cosd(lat1).*cosd(lon1);cosd(lat1).*sind(lon1);sind(lat1)];

%Görlitz

% Geozentrische Satellitenposition zu den Zeiten t1, t2
x = SatData1.xyz(1,:);
y = SatData1.xyz(2,:);
z = SatData1.xyz(3,:);

% Geozentrische Ortsvektoren
vecx   = [x;y;z];
r1     = vecnorm(vecx);

% Koordinatentransformation
kend = length(theta1);
for k= 1:kend
    xS(:,k) = mtimes(R_z(theta1(k)),vecx(:,k));
end
XS       = xS - vecXB; 
distance = vecnorm(XS);
S        = mtimes(U,XS);
 
Az       = atan2d(S(1,:),S(2,:));
hoehe    = atan2d(S(3,:),(sqrt(S(1,:).^2+S(2,:).^2)));

for k = 1:length(Az)
   if hoehe(k)< 0 
       Az(k) = NaN;
       distance1(k) = NaN;
   end
end

titlestr = strcat(' Topozentr. Beobachtung (',LocStr(1,:),')');
figure('name',titlestr);
subplot(2,1,1)
yyaxis left
plot(datetime(T_vector1,'ConvertFrom','juliandate'),hoehe,...
    'LineWidth',2);
hold on
ylabel('Hoehe in °')
ylim([0,90]);
yyaxis right
ylabel('Azimut in °')
p(1)=plot(datetime(SatData1.Time,'ConvertFrom','juliandate'),Az,...
    'LineWidth',2);
hold on
grid on;
hp1 = title(strcat('Höhe und Azimut (',LocStr(1,:),')'), 'FontSize',12);
set(hp1,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',14);

subplot(2,1,2)
plot(datetime(SatData1.Time,'ConvertFrom','juliandate'),...
     distance,'color',Colors(3,:),'LineWidth',2);
hold on
ylabel('Entfernung in km')
% ylim([0,1.1*max(max(distance1),max(distance2))]);
grid on;
hp1 = title(strcat('Topoz. Entfernung (',LocStr(1,:),')'), 'FontSize',12);
set(hp1,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',14);

%% Ende Programm


