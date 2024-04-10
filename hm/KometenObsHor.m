% -------------------------------------------------------------------------
% KometenObsHor.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die Kometenbahnen als Funktion der Zeit in Erdnaehe und die 
% Beobachtungsdaten (Azimut und Hoehe) ueber einen Zeitraum von mehreren 
% Tagen.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
% Auswahl des Kometen
NrK=2;  % Neowise

% Koordinaten Oberkochen (Beobachtungsort) 
lambda=deg2rad(10.0990193);
phi=deg2rad(48.7854368);

GM      = 2.95479E-04;          % µG [AE^3/d^2]
c_l     = 173.14;               % Lichtgeschwindigkeit [AE/d]
DT      = 70/86400;             % DeltaT in Tagen;

%--------------------------------------------------------------------------
% Beginn Rechnung 1
%% 
% Zeitraum fuer Berechnung Alt-AZ ueber mehrere Tage (nach Mitternacht)

dt1 = datetime('2020-07-01 03:00:00');
dt2 = datetime('2020-07-31 03:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
T2 = juliandate(dt2); % Julianisches Datum  
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');
TimeStr1 = strjoin(["Zeit : ",string(dt1,'HH:mm'), " (UT)"]);

Aequi = 'Datum';
% Einlesen der Bahnparameter der Planeten und deren Ableitungen 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);
epsErad = deg2rad(EpsErde(T1)); % Ekliptikneigung 2020


Nges            = round((T2-T1))+1;     % Anzahl betrachtete Tage  
T               = linspace(T1,T2,Nges);
MJuDa_vector    = linspace(MJuDa1,MJuDa2,Nges);

% Einstellungen Schrittweiten 
StepsperDay=1;
NrStepsMarker=StepsperDay*2;
NrStepsLabel=NrStepsMarker*2;
MaxStepsLabel=floor(Nges/NrStepsLabel);
 
% Einlesen der Bahnparameter der Kometen 
BaPaK=KometParameter('Kometen.csv');

% Berechnung der heliozentrischen Koordinaten der Erde und des Kometen 
% nach Keplerloesung von T1 bis T2 

Erde   = PlanetPQR(T, BaPa, BaPadot, 3);
Komet  = KometPQR(GM, T, BaPaK,NrK);

OmegaPn=deg2rad(BaPaK.OmegaP(NrK));
omegaPn=deg2rad(BaPaK.omegaP(NrK));
iPn=deg2rad(BaPaK.iP(NrK));

% Berechnung Abstand Erde-Komet
rek=Komet.xyz-Erde.xyz;
rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));

% Beruecksichtigung Lichtlaufzeit
rek_korr= rek-rekbetr.*Komet.v/c_l;

% Koordinatentransformation in aequatoriale Koordinaten
Komet.equ = CalcAnglesfromXYZ(mtimes(R_x(-epsErad),rek_korr));
[RAS,DECS]=KeplerSonne(T,epsErad); 

% aequatoriale und horizontale Koordianten
RA(:)=rad2deg(wrapToPi(Komet.equ(2,:)));
DEC(:)=rad2deg(Komet.equ(3,:));
RAS(:)=rad2deg(RAS);
DECS(:)=rad2deg(DECS);

% Berechnung des Stundenwinkels und Alt-Az
Komet.Tau = GMST(MJuDa_vector)+lambda-deg2rad(RA);
SonneTau=GMST(MJuDa_vector)+lambda-deg2rad(RAS);
Alt= asind(sin(phi).*sind(DEC)+cos(phi)*cosd(DEC).*cos(Komet.Tau));
Az = wrapTo360(atan2d(cosd(DEC).*sin(Komet.Tau), ...
    (cosd(DEC).*cos(Komet.Tau).*sin(phi)-sind(DEC)*cos(phi))))-180;

AltS= asind(sin(phi).*sind(DECS)+cos(phi)*cosd(DECS).*cos(SonneTau));
AzS = wrapTo360(atan2d(cosd(DECS).*sin(SonneTau), ...
    (cosd(DECS).*cos(SonneTau).*sin(phi)-sind(DECS)*cos(phi))))-180;

%%
%--------------------------------------------------------------------------
% Beginn Rechnung 2
% Zeitraum fuer Berechnung Alt-AZ ueber mehrere Tage (vor Mitternacht)

dt3 = datetime('2020-07-01 20:00:00');
dt4 = datetime('2020-07-31 20:00:00');

T3 = juliandate(dt3); % Julianisches Datum
MJuDa3 = juliandate(dt3,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr3 =  string(dt3,'dd.MM.yyyy');
T4 = juliandate(dt4); % Julianisches Datum  % Bedeckungsdatum
MJuDa4 = juliandate(dt4,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr4 =  string(dt4,'dd.MM.yyyy');
TimeStr3 = strjoin(["Zeit : ",string(dt3,'HH:mm'), " (UT)"]);

Nges=round((T4-T3))+1;
T=linspace(T3,T4,Nges);
MJuDa_vector2=linspace(MJuDa3,MJuDa4,Nges);
  
% Berechnung der heliozentrischen Koordinaten der Erde und des Kometen 
% nach Keplerloesung von T3 bis T4 

Erde   = PlanetPQR(T, BaPa, BaPadot, 3);
Komet2 = KometPQR(GM, T, BaPaK,NrK);

% Berechnung Abstand Erde-Komet
rek=Komet2.xyz-Erde.xyz;
rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));
[minrekbetr,ErdNaehe] = min(rekbetr);

% Beruecksichtigung Lichtlaufzeit
rek_korr= rek-rekbetr.*Komet2.v/c_l;

% Koordinatentransformation in aequatoriale Koordinaten
Komet2.equ = CalcAnglesfromXYZ(mtimes(R_x(-epsErad),rek_korr)); 
[RAS,DECS]=KeplerSonne(T,epsErad); 

% aequatoriale und Horizontale Koordianten
RA(:)=rad2deg(wrapToPi(Komet2.equ(2,:)));
DEC(:)=rad2deg(Komet2.equ(3,:));
RAS(:)=rad2deg(RAS);
DECS(:)=rad2deg(DECS);

% Berechnung des Stundenwinkels und Alt-Az
Komet2.Tau = GMST(MJuDa_vector2)+lambda-deg2rad(RA);
SonneTau=GMST(MJuDa_vector2)+lambda-deg2rad(RAS);
Alt2= asind(sin(phi).*sind(DEC)+cos(phi)*cosd(DEC).*cos(Komet2.Tau));
Az2 = wrapTo360(atan2d(cosd(DEC).*sin(Komet2.Tau), ...
    (cosd(DEC).*cos(Komet2.Tau).*sin(phi)-sind(DEC)*cos(phi))))-180;

AltS2= asind(sin(phi).*sind(DECS)+cos(phi)*cosd(DECS).*cos(SonneTau));
AzS2 = wrapTo360(atan2d(cosd(DECS).*sin(SonneTau), ...
    (cosd(DECS).*cos(SonneTau).*sin(phi)-sind(DECS)*cos(phi))))-180;

%--------------------------------------------------------------------------
%%
%Graphische Ausgabe

% Darstellung horizontaler Koordinaten
header2=' Bahnen von Kometen';
figure('Name',header2);
header2 = cellstr(strcat(BaPaK.Name(NrK)));

% Plot Bahn des Kometen 
plot(Az,Alt,':+','MarkerIndices',1:NrStepsMarker:length ...
    (Komet.xyz(1,:)),'LineWidth',2,'Color',Colors(NrK,:));
hold on;
plot(Az2,Alt2,':+','MarkerIndices',1:NrStepsMarker:length ...
    (Komet.xyz(1,:)),'LineWidth',2,'Color',Colors(NrK,:));
% Plot Bahn der Sonne 
plot(AzS,AltS,':+','MarkerIndices',1:NrStepsMarker*2:length ...
    (Komet.xyz(1,:)),'LineWidth',2,'Color',Colors(10,:));
plot(AzS2,AltS2,':+','MarkerIndices',1:NrStepsMarker*2:length ...
    (Komet.xyz(1,:)),'LineWidth',2,'Color',Colors(10,:));

% Einstellungen der Achsen 
xlim([-90 60]);
ylim([-20 40]);
hold on;

% Markierung der Nulllinie 
x=linspace(-150,150,10);
y=zeros(10);
plot(x,y,'-.','LineWidth',2,'Color',Colors(5,:));

% Beschriftung der Zeitpunkte 
for k=1:MaxStepsLabel+1
       DateObs = datetime(Komet.Time(1+(k-1)*NrStepsLabel), ...
           'convertfrom','juliandate');
       mylabels(k,:) = string(DateObs,'dd-MM-yyyy');
       Xlab(k)=Az((k-1)*NrStepsLabel+1);
       Ylab(k)=Alt((k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,Ylab,mylabels,'S',0.02,0,'FontSize',8, ...
    'Color',Colors(NrK,:));
LabelPoints(Az(Nges),Alt(Nges),TimeStr1,'S',0.05,0,'FontSize',12, ...
    'Color',Colors(NrK,:));
for k=1:MaxStepsLabel+1
       DateObs = datetime(Komet2.Time(1+(k-1)*NrStepsLabel), ...
           'convertfrom','juliandate');
       mylabels(k,:) = string(DateObs,'dd-MM-yyyy');
       Xlab(k)=Az2((k-1)*NrStepsLabel+1);
       Ylab(k)=Alt2((k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,Ylab,mylabels,'S',0.02,0,'FontSize',8, ...
    'Color',Colors(NrK,:));
LabelPoints(Az2(Nges),Alt2(Nges),TimeStr3,'S',0.05,0,'FontSize',12, ...
    'Color',Colors(NrK,:));
LabelPoints(AzS2(Nges),AltS2(Nges),'Sonne','W',0.05,0,'FontSize',12, ...
    'Color',Colors(10,:));
LabelPoints(AzS(Nges),AltS(Nges),'Sonne','W',0.05,0,'FontSize',12, ...
    'Color',Colors(10,:));

grid on;
grid minor;
txt=strjoin([header2, ' ', DatumStr1, " - ", DatumStr2]);
title(txt,'FontSize',12);
ylabel('Hoehe in °');
xlabel('Azimut in °');
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
