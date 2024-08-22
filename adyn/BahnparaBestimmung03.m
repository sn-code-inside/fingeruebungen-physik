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
%  
% Bestimmung der Bahnparameter aus topozentrischen 
% Beobachtungsdaten über das Gauß-Verfahren.
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

% Einheitsvektoren, Transformationsmatrix Topozentrisch<-->Geozentrisch KOS
eO = [-sin(lambdaB);cos(lambdaB);0];
eN = [-sin(phiB).*cos(lambdaB);-sin(phiB).*sin(lambdaB);cos(phiB)];
eZ = [cos(phiB).*cos(lambdaB);cos(phiB).*sin(lambdaB);sin(phiB)];
U  = [eO,eN,eZ]';
UT = transpose(U);


% Beobachtungszeiten
dt1 = datetime(['2020-01-01 00:00:00']);
T1  = juliandate(dt1); % Julianisches Datum  ET
dt2 = datetime('2020-01-01 00:15:00');
T2  = juliandate(dt2); % Julianisches Datum  ET
dt3 = datetime('2020-01-01 00:30:00');
T3  = juliandate(dt3); % Julianisches Datum  ET
NPoints1     = 151;
NPoints2    =  301;
T_vector1   = linspace(T1,T2,NPoints1);
T_vector2   = linspace(T1,T3,NPoints2);

%% Bahnparameter Satellit(Vorwärtsrechnung)

aP      = RE+2000;
eP      = 0.1;
iP      = 95; 
OmegaP  = 135;
omegaP  = 10;
TPn     = 2*pi*sqrt(aP^3/GME)/86400;
MP      = 10;

SatHEO.BaPa = [aP , eP, iP , MP, OmegaP, omegaP];
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
lon2d(:)    = rad2deg(lon2);


fprintf('\n');
fprintf('Bahnparameter für Bahnberechnung nach Kepler\n');
fprintf('\n Zeit  = %s', dt1);
fprintf('\n');
fprintf('\n a     = %8.2f km  (Große Halbachse)', aP);
fprintf('\n e     = %8.5f     (Exzentrizität)', eP);
fprintf('\n i     = %8.2f °   (Inklination)', iP);
fprintf('\n Omega = %8.2f °   (Rektasz. aufst. Knoten)', OmegaP);
fprintf('\n omega = %8.2f °   (Argument Perigäum)', omegaP);
fprintf('\n M     = %8.2f °   (Mittlere Anomalie Beobachtungsstart)', MP);
fprintf('\n TP    = %8.2f min (Umlaufzeit)', TPn*24*60);
fprintf('\n TP    = %8.2f h   (Umlaufzeit)', TPn*24);
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
geolimits([-90 90],[-180 +90])
title(titlestr);
set(gca,'FontSize',14);


%% Berechnung topozentrischer Beobachtungsdaten

% Geozentrische Satellitenposition zu den Zeiten t1, t2
MP1= 61;
MP2=141;
xA = SatData1.xyz(1,MP1);
yA = SatData1.xyz(2,MP1);
zA = SatData1.xyz(3,MP1);
xE = SatData1.xyz(1,MP2);
yE = SatData1.xyz(2,MP2);
zE = SatData1.xyz(3,MP2);

T1m  = SatData1.Time(MP1);
T2m  = SatData1.Time(MP2);
dtm1 = datetime(T1m,'ConvertFrom','juliandate') ;
dtm2 = datetime(T2m,'ConvertFrom','juliandate') ;

% Geozentrische Ortsvektoren
vecx1  = [xA;yA;zA];
vecx2  = [xE;yE;zE];
r1m     = norm(vecx1);
r2m     = norm(vecx2);

% Koordinatentransformation
xS1   = R_z(theta1(MP1))*vecx1;
xS2   = R_z(theta1(MP2))*vecx2;
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
fprintf('\n Topozentrische Beobachtung: ')
fprintf('\n Länge            = %+6.2f°',  lambdaBd);
fprintf('\n Breite           = %+6.2f°',  phiBd);
fprintf('\n ');
fprintf('\n Beobachtung für %s UTC \n', dtm1);
fprintf('\n Kart. Koordinaten  x1  = [%+10.3f  %+10.3f  %+10.3f] km ', xA, yA, zA);
fprintf('\n Azimut           = %+6.2f°', Az1);
fprintf('\n Höhe             = %+6.2f° ', hoehe1);
fprintf('\n Abstand          = %+6.2f  km', distance1);
fprintf('\n ');
fprintf('\n Beobachtung für %s UTC \n', dtm2);
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
% r1     = vecnorm(vecx);

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
yticks([0,30,60,90]);

yyaxis right
ylabel('Azimut in °')
p(1)=plot(datetime(SatData1.Time,'ConvertFrom','juliandate'),Az,...
    'LineWidth',2);
hold on
grid on;
ylim([0,180]);
yticks([0 45 90 135 180]);
hp1 = title(strcat('Höhe und Azimut (',LocStr(1,:),')'), 'FontSize',12);
set(hp1,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',14);

subplot(2,1,2)
plot(datetime(SatData1.Time,'ConvertFrom','juliandate'),...
     distance,'color',Colors(3,:),'LineWidth',2);
hold on
ylabel('Entfernung in km')
ylim([0,aP]);
grid on;
ylim([0 6000])
yticks([0,2000,4000,6000])
hp1 = title(strcat('Topoz. Entfernung (',LocStr(1,:),')'), 'FontSize',12);
set(hp1,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',14);


%% Rückwärtsrechnung nach Gauss

dist1m    = distance(MP1);
dist2m    = distance(MP2);

% Berechnung der kartesischen, geozentrischen Koordinaten aus 
% topozentrischen Beobachtungsdaten zu T1m, T2m

[RO, RN, RZ] = sph2cart(deg2rad(90-Az1),deg2rad(hoehe1),dist1m); 
%Achtung x-Achse = Ost, 90° geg AZ gedreht
vecS1m       = [RO; RN; RZ];
thet1md      = GMSTsat(T1m);
thetm1       = deg2rad(thet1md);
vecx1m       = mtimes(transpose(R_z(thetm1)),(mtimes(UT,vecS1m)+vecXB));

[RO, RN, RZ] = sph2cart(deg2rad(90-Az2),deg2rad(hoehe2),dist2m);
vecS2m       = [RO; RN; RZ];
thet2md      = GMSTsat(T2m);
thetm2       = deg2rad(thet2md);
vecx2m       = mtimes(transpose(R_z(thetm2)),(mtimes(UT,vecS2m)+vecXB));


% % Oder direkte Benutzung geozentrische Satellitenposition zu Zeiten T1m,T2m
vecx1m = [xA;yA;zA];
vecx2m = [xE;yE;zE];


%--------------------------------
r1m     = norm(vecx1m);
r2m     = norm(vecx2m);
vece1 = vecx1m/r1m;
vecr0 = vecx2m-dot(vecx2m,vece1)*vece1;
vece0 = vecr0/norm(vecr0);
r0    = norm(vecr0);
W     = cross(vece1,vece0);

inc   = atan2d(sqrt(W(1)^2+W(2)^2),W(3));
Omega = wrapTo360(atan2d(W(1),-W(2)));
u1    = atan2d(vecx1(3),(-vecx1(1)*W(2)+vecx1(2)*W(1)));

tau   = sqrt(GME)*(T2m-T1m)*86400;
NN    = sqrt(2*(r1m*r2m+dot(vecx1,vecx2)));
mpara = tau^2/NN^3;
kpara = (r1m+r2m)/2/NN -0.5;

% ------------------------------------------------------------------------
% Ansatz
% Nullstelle der nichtlinearen Funktion funF(eta) in eta mit 
% der Routine fzero

% [RO, RN, RZ] = sph2cart(pi/2-Az1,hoehe1,r1); %Achtung x-Achse = Ost, 90° geg AZ
% vecS1m         = [RO; RN; RZ];
% thet1md        = GMSTsat(T1);
% thet1m         = deg2rad(thet1md);
% vecx1M         = mtimes(transpose(R_z(thet1m)),(mtimes(UT,vecS1m)+vecXB));
% r1             = norm(vecx1M);
% [RO, RN, RZ]   = sph2cart(pi/2-Az2,hoehe2,r2);
% vecS1m         = [RO; RN; RZ];
% thet2md        = GMSTsat(T2);
% thet2m         = deg2rad(thet2md);
% vecx2M         = mtimes(transpose(R_z(thet2m)),(mtimes(UT,vecS1m)+vecXB));
% r2             = norm(vecx2M);
% 
% vece1 = vecx1M/r1;
% vecr0 = vecx2M-dot(vecx2M,vece1)*vece1;
% vece0 = vecr0/norm(vecr0);
% r0    = norm(vecr0);
% W     = cross(vece1,vece0);

% inc   = atan2d(sqrt(W(1)^2+W(2)^2),W(3));
% Omega = wrapTo360(atan2d(W(1),-W(2)));
% u1    = atan2d(vecx1M(3),(-vecx1M(1)*W(2)+vecx1M(2)*W(1)));
% 
% tau   = sqrt(GME)*(T2m-T1m)*86400;
% NN    = sqrt(2*(r1*r2+dot(vecx1M,vecx2M)));
% mpara = tau^2/NN^3;
% kpara = (r1+r2)/2/NN -0.5;

% ------------------------------------------------------------------------
% Ansatz  Funktion funF(eta) in eta mit 
% der Routine fzero

Dups  = acosd(dot(vecx1m,vecx2m)/r1m/r2m);
myfunF = @(eta,P1) 1- eta + P1.m./eta^2*(4*asin(sqrt(P1.m/eta^2-P1.k))- ...
         sin(4*asin(sqrt(P1.m/eta^2-P1.k))))/...
         sin(2*asin(sqrt(P1.m/eta^2-P1.k)))^3;
P1.m = mpara;
P1.k = kpara;
funF = @(x) myfunF(x,P1);  % function of x alone
etaF = fzero(funF,-1000);

p     = (etaF*r1m*r2m*sind(Dups)/tau)^2;
ups1  = atan2d(((p/r1m-1)*dot(vecx2m,vece1)/r2m-(p/r2m-1)),r0*(p/r1m-1)/r2m);
exz   = (p/r1m-1)/cosd(ups1);
if exz > 1
    disp('   ')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('Numerische Lösung ergibt keine Ellipse')
    disp('Parameter zu ungenau oder Zeit/Winkelabstand zu gering')
    return
end
a     = p/(1-exz^2);
omega = wrapTo360(u1-ups1);
ups2  = acosd((p/r2m-1)/exz);
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
fprintf('\n M     = %9.3f °   (Mittlere Anomalie Messung 1)', M);
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
 
eta0  = 12/22+10*sqrt(1+(44/9)*mpara/(kpara+5/6))/22;
eta1  = eta0+0.1;
% Erste Iteration
eta2  = funciter(eta1, eta0, mpara,kpara);
% Zweite Iteration
eta0  = eta1;
eta1  = eta2;
eta2  = funciter(eta1, eta0, mpara,kpara);
% Dritte Iteration
eta0  = eta1;
eta1  = eta2;
eta2  = funciter(eta1, eta0, mpara,kpara);

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
