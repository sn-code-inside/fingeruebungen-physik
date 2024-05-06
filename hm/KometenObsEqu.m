% -------------------------------------------------------------------------
% KometenObsEqu.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die Kometenbahnen in Erdnaehe als Funktion der Zeit und die 
% Beobachtungsdaten (Rektaszension und Deklination) unter Benutzung von 
% parabolischen Naeherungen.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
% Auswahl des Kometen
NrK=2;

% Initialisierung
% Einlesen der Bahnparameter der Kometen 
BaPaK=KometParameter('Kometen.csv');
T_Per   = BaPaK.T0(NrK);
GM      = 2.95479E-04;     % µG [AE^3/d^2]
c_l     = 173.14;          % Lichtgeschwindigkeit [AE/d]

% Zeitraum 
dt1 = datetime('2020-07-01 00:00:01');
dt2 = datetime('2020-08-01 00:00:01');

T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
T2 = juliandate(dt2); % Julianisches Datum  
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');
TimeStr = string(dt1,'HH:mm') + " ET";
eps0 = deg2rad(EpsErde(T1));

Aequi = 'Datum';
% Einlesen der Bahnparameter der Planeten und deren Ableitung 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);

% Einstellung der Schrittweiten 
Nges=floor(T2-T1)*24;
T=linspace(T1,T2,Nges);
StepsperDay=24;
NrStepsMarker=StepsperDay;
NrStepsLabel=NrStepsMarker*5;
MaxStepsLabel=floor(Nges/NrStepsLabel);
   

%--------------------------------------------------------------------------
%%
% Beginn Rechnung

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
[minrekbetr,ErdNaehe] = min(rekbetr);

% Beruecksichtigung Lichtlaufzeit
rek_korr= rek-rekbetr.*Komet.v/c_l;

% Koordinatentransformation in aequatoriale Koordinaten 
Komet.equ = CalcAnglesfromXYZ(mtimes(R_x(-eps0),rek_korr)); 
[RAS,DECS]=KeplerSonne(T,eps0); 

% aequatoriale und horizontale Koordinaten
RA(:)=rad2deg(wrapToPi(Komet.equ(2,:)));
DEC(:)=rad2deg(Komet.equ(3,:));
RAS(:)=rad2deg(RAS);
DECS(:)=rad2deg(DECS);

%--------------------------------------------------------------------------
%%
% Graphische Ausgabe

% Darstellung aequatorialer Koordinaten
header2='Kometenbahn Dec-RA';
figure('Name',header2);
plot(RA,DEC,'Color', Colors(NrK,:));
hold on;
plot(RAS,DECS,'Color', Colors(NrK,:));
plot(RA,DEC,':+','MarkerIndices',1:NrStepsMarker:length ...
    (Komet.xyz(1,:)),'LineWidth',1,'Color',Colors(NrK,:));
plot(RAS,DECS,':+','MarkerIndices',1:NrStepsMarker*5:length ...
    (Komet.xyz(1,:)),'LineWidth',1,'Color',Colors(4,:));
hold on;
% Beschriftung der Zeitpunkte 
for k=1:MaxStepsLabel+1
       DateObs = datetime(Komet.Time(1+(k-1)*NrStepsLabel), ...
           'convertfrom','juliandate');
       mylabels(k,:) = string(DateObs,'dd-MM-yyyy');
       Xlab(k)=RA((k-1)*NrStepsLabel+1);
       Ylab(k)=DEC((k-1)*NrStepsLabel+1);
       XlabS(k)=RAS((k-1)*NrStepsLabel+1);
       YlabS(k)=DECS((k-1)*NrStepsLabel+1);
end

% Einfuegen der Beschriftungen innerhalb des Schaubilds 
LabelPoints(RA(Nges),DEC(Nges),header2,'SW',0.04,0,'FontSize',12, ...
    'Color',Colors(NrK,:));
LabelPoints(Xlab,Ylab,mylabels,'E',0.02,0,'FontSize',8, ...
    'Color',Colors(NrK,:));
LabelPoints(XlabS,YlabS,mylabels,'N',0.01,0,'FontSize',8, ...
    'Color',Colors(4,:));
LabelPoints(RAS(Nges),DECS(Nges),'Sonne','SW',0.04,0,'FontSize',12, ...
    'Color',Colors(4,:));

grid on;
grid minor;

% Titel 
header2 = cellstr(strcat(BaPaK.Name(NrK)));
header3 = sprintf("%s %s ",' - taeglich um',TimeStr);
ttl=title(strcat(header2, header3));
ttl.FontSize=14;
ylabel('Dec in °');
xlabel('RA in °');
% Einstellungen der Achsen 
xlim([80 200]);
ylim([10 50]);
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
