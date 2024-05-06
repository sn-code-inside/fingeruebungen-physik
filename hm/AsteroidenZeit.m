% -------------------------------------------------------------------------
% AsteroidenZeit.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahn eines Asteroiden als Funktion der Zeit 
% unter fallweiser Benutzung einer parabolischen Naeherung.
% (Aequinoktium und Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Parameter einlesen
T1 = juliandate(datetime(2000,1,1,12,0,0));
Aequi = 'J2000';
% Einlesen der Bahnparameter der Planeten und deren Ableitung
[BaPa,BaPadot]=OrbitParameter(T1, Aequi); 
% Einlesen der Bahnparameter der Kometen und Asteroiden 
BaPaK=KometParameter('Kometen.csv'); 

GM      = 2.95479E-04;     % µG [AE^3/d^2]
c_l     = 173.14;          % Lichtgeschwindigkeit [AE/d]

% Auswahl des Asteroiden (Nr 9, 10 oder 11 eingeben)
AK=10;
AE=2;  
% Zeitraum und Darstellung Perihelnaehe
dt1 = datetime('2021-03-01 00:00:00');
dt2 = datetime('2021-06-01 00:00:00');
T1 = juliandate(dt1);
T2 = juliandate(dt2);

Nges=10000; % Zahl der Stuetzstellen im Intervall T1 - T2
StepsperDay=round(Nges/(T2-T1)); % Zahl der Stuetzstellen pro Tag

% Input markierte Tage
StepsMarker= 5; %Zahl der markierten Tage 
% (Richtwert 1 bei T2-T1=1 Monat,Richtwert 24 bei T2-T1=1 Jahr usw.)
NrStepsMarker=StepsMarker*StepsperDay; %Zahl der markierten Tage

% Input angezeigte Tage 
StepsLabel=2; %Zahl der angezeigten Tage
%((Richtwert 2 bei T2-T1=1 Monat,Richtwert 5 bei T2-T1=1 Jahr usw.)
NrStepsLabel=StepsLabel*NrStepsMarker; % Zahl der angezeigten Tage
MaxStepsLabel=floor(Nges/NrStepsLabel)+1;
T=linspace(T1,T2,Nges); % Berechnungsintervall

DatumStr1 =  string(dt1,'dd.MM.yy');
DatumStr2 =  string(dt2,'dd.MM.yy');

T_Per     = BaPaK.T0(AK);           % Perihelzeit des Asteroiden
Ax_L      = 1.5*BaPaK.qP(AK);       % Achsenlimits auf 2*Perihelabstand


%%
%--------------------------------------------------------------------------
% Beginn Rechnung

% Berechnung der Erdposition zum Perihel des Asteroiden
ErdePeri=PlanetPQR(T_Per, BaPa, BaPadot, 3);

% Berechnung der heliozentrischen Koordinaten der Erde und des Mars
% von T1 bis T2 

for k=3:4
    PlanetPos(k)=PlanetPQR(T, BaPa, BaPadot, k);
end

% Berechnung der heliozentrischen Koordinaten 
% des Asteroiden nach Loesung mit Stumppfschen Funktionen
% von T1 bis T2 

Komet  = KometPQR(GM, T, BaPaK, AK);

OmegaPn=deg2rad(BaPaK.OmegaP(AK));
omegaPn=deg2rad(BaPaK.omegaP(AK));
iPn=deg2rad(BaPaK.iP(AK));

% PQR Multiplikation Perihel
yin=[BaPaK.qP(AK);0;0];
Perihel.xyz=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);

% Berechnung der Planetenbahnen Erde und Mars Ellipsengleichung
Planets = PlanetOrbits(BaPa);

% Berechnung Abstand Erde-Komet
rek=Komet.xyz-PlanetPos(3).xyz;
rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));
[minrekbetr,ErdNaehe] = min(rekbetr);

% Beruecksichtigung Lichtlaufzeit
rek_korr= rek-rekbetr.*Komet.v/c_l;

% aequatoriale Koordinaten
T0=(T1-2451545)/36525;
eps0=deg2rad(23.43929111);
% Koordinatentrafo in Equ-System mit Beruecksichtigung Lichtlaufzeit
Komet.equ = CalcAnglesfromXYZ(mtimes(R_x(-eps0),rek_korr)); 
% Koordinatentrafo in Equ-System ohne Beruecksichtigung Lichtlaufzeit
% Komet.equ = CalcAnglesfromXYZ(mtimes(R_x(-eps0),rek)); 

%--------------------------------------------------------------------------
%%
% Graphische Ausgabe

%__________________________________________________________________________
%
% Figure 1: Ekliptikebenen - Darstellung
%__________________________________________________________________________

header1='Asteroidenbahnen';
figure('Name',header1);

% x-y-Ebene
subplot(1,2,1);         % linkes Schaubild in Figure 1 
% Plot Asteroid
plot(Komet.xyz(1,:),Komet.xyz(2,:),':+','MarkerIndices', ...
    1:NrStepsMarker:length(Komet.xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(AK,:));
hold on;
for k=1:MaxStepsLabel     
    dt4 = datetime(Komet.Time((k-1)*NrStepsLabel+1),'ConvertFrom', ...
        'juliandate');
    mylabels(k,:)= string(dt4,'dd-MM-yy');
end
for k=1:MaxStepsLabel
    Xlab(k)=Komet.xyz(1,(k-1)*NrStepsLabel+1);
    Ylab(k)=Komet.xyz(2,(k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,Ylab,mylabels,'E',0.01,1,'FontSize',8,'Color', ...
    Colors(AK,:));
 
% Sonne, Perihel , Richtung Fruehlingspunkt
SonneFP;
x=linspace(0,Perihel.xyz(1),10);
y=linspace(0,Perihel.xyz(2),10);
p(10)=plot(x,y,'Color',Colors(10,:));
p(10).LineWidth=2;

% Position Erdnaehe
x=linspace(PlanetPos(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe),10);
y=linspace(PlanetPos(3).xyz(2,ErdNaehe),Komet.xyz(2,ErdNaehe),10);
p(10)=plot(x,y,'Color',Colors(3,:));
p(10).LineWidth=2;
plot(PlanetPos(3).xyz(1,ErdNaehe),PlanetPos(3).xyz(2,ErdNaehe),'o', ...
    'Color', Colors(AK,:),'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
h=LabelPoints(PlanetPos(3).xyz(1,ErdNaehe), ...
    PlanetPos(3).xyz(2,ErdNaehe),'Erdnaehe','E',0.02,'FontSize',8, ...
    'Color',Colors(3,:));

% Planetenbahnen
for iPlot = 3:4
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
        Colors(iPlot,:));
end

% Position Planeten an Anfang und Ende des Intervalls
plot(ErdePeri.xyz(1,1),ErdePeri.xyz(2,1),'+','Color', Colors(3,:), ...
    'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
h=LabelPoints(ErdePeri.xyz(1,1), ErdePeri.xyz(2,1), ...
    'Erde@Asteroid-Perihel','S',0.02,'FontSize',8,'Color',Colors(3,:));
for iPlot = 3:4
    plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1),'+', ...
        'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1), ...
        Planets(iPlot).Name+'@T1','E',0.02,'FontSize',8,'Color', ...
        Colors(iPlot,:));
    plot(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(2,Nges),'+', ...
        'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,Nges), ...
        PlanetPos(iPlot).xyz(2,Nges),Planets(iPlot).Name+'@T2','E', ...
        0.02,'FontSize',8,'Color',Colors(iPlot,:));
end

header2 = strcat(BaPaK.Name(AK));
title(header2);
% Einstellungen Achsen 
xlim([-Ax_L Ax_L]);
ylim([-Ax_L Ax_L]);
axis square;
grid on;
grid minor;
xlabel('x in AE');
ylabel('y in AE');

% x-z-Ebene

subplot(1,2,2);         % rechtes Schaubild in Figure 1 
hold on;
% Plot Asteroid
plot(Komet.xyz(1,:),Komet.xyz(3,:),'-+','MarkerIndices', ...
    1:NrStepsMarker:length(Komet.xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(AK,:));
for k=1:MaxStepsLabel     
    dt4 = datetime(Komet.Time((k-1)*NrStepsLabel+1),'ConvertFrom', ...
        'juliandate');
    mylabels(k,:)= string(dt4,'dd-MM-yy');
end
for k=1:MaxStepsLabel
    Xlab(k)=Komet.xyz(1,(k-1)*NrStepsLabel+1);
    Ylab(k)=Komet.xyz(3,(k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,Ylab,mylabels,'E',0.01,1,'FontSize',8,'Color', ...
    Colors(AK,:));

% Sonne, Perihel , Richtung Fruehlingspunkt
SonneFP;
x=linspace(0,Perihel.xyz(1),10);
y=linspace(0,Perihel.xyz(3),10);
p(10)=plot(x,y,'Color',Colors(10,:));
p(10).LineWidth=2;

% Position Erdnaehe
x=linspace(PlanetPos(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe),10);
z=linspace(PlanetPos(3).xyz(3,ErdNaehe),Komet.xyz(3,ErdNaehe),10);
p(10)=plot(x,z,'Color',Colors(3,:));
p(10).LineWidth=2;
plot(PlanetPos(3).xyz(1,ErdNaehe),PlanetPos(3).xyz(3,ErdNaehe),'o', ...
    'Color', Colors(AK,:),'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
h=LabelPoints(PlanetPos(3).xyz(1,ErdNaehe), ...
    PlanetPos(3).xyz(3,ErdNaehe),'Erdnaehe','E',0.02,'FontSize',8, ...
    'Color',Colors(3,:));

% Planetenbahnen
for iPlot = 3:4
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(3,:),'Color', ...
        Colors(iPlot,:));
end
% Position Planeten Anfang und Ende Intervall, Erde am Perihel
plot(ErdePeri.xyz(1,1),ErdePeri.xyz(3,1),'+','Color', Colors(3,:), ...
    'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
h=LabelPoints(ErdePeri.xyz(1,1), ErdePeri.xyz(3,1), ...
    'Erde@Asteroid-Perihel','S',0.02,'FontSize',8,'Color',Colors(3,:));
for iPlot = 3:4
    plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(3,1),'+', ...
        'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,1), ...
        PlanetPos(iPlot).xyz(3,1),Planets(iPlot).Name+'@T1','E',0.02, ...
        'FontSize',8,'Color',Colors(iPlot,:));
    plot(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(3,Nges), ...
        '+','Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,Nges), ...
        PlanetPos(iPlot).xyz(3,Nges),Planets(iPlot).Name+'@T2','E', ...
        0.02,'FontSize',8,'Color',Colors(iPlot,:));
end

header2 = strcat(BaPaK.Name(AK));
title(header2);

% Einstellungen Achsen 
xlim([-Ax_L Ax_L]);
ylim([-Ax_L Ax_L]);
axis square;
grid on;
grid minor;
xlabel('x in AE');
ylabel('z in AE');

%__________________________________________________________________________
%
% Figure 2: 3-Dimensionale Darstellung
%__________________________________________________________________________


figure('Name',header1);
p=plot3(Komet.xyz(1,:),Komet.xyz(2,:),Komet.xyz(3,:),'Color',Colors(AK,:));
p.LineWidth=2;
hold on
axis equal

%Sonne und Richtung zum Fruehlingspunkt
SonneFP;
% Iteration durch die Planeten (hier: Erde und Mars) fuer Planetenbahnen
for iPlot= 3:4
  plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:), ...
      Planets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
  hold on;
  axis square;
end

% Position Planeten Anfang und Ende Intervall, Erde am Perihel
plot3(ErdePeri.xyz(1,1),ErdePeri.xyz(2,1),ErdePeri.xyz(3,1),'o', ...
    'Color', Colors(3,:),'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
text(ErdePeri.xyz(1,1), ErdePeri.xyz(2,1), ErdePeri.xyz(3,1), ...
    'Erde@Asteroid-Perihel','VerticalAlignment','top','FontSize',8, ...
    'Color',Colors(3,:));

% Iteration durch die Planeten (hier: Erde und Mars) fuer Positionen 
for iPlot = 3:4
    plot3(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1), ...
        PlanetPos(iPlot).xyz(3,1),'+','Color', Colors(iPlot,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
    text(PlanetPos(iPlot).xyz(1,1), PlanetPos(iPlot).xyz(2,1), ...
        PlanetPos(iPlot).xyz(3,1),Planets(iPlot).Name+'@T1', ...
        'VerticalAlignment','bottom','FontSize',8,'Color',Colors(iPlot,:));
    plot3(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(2,Nges), ...
        PlanetPos(iPlot).xyz(3,Nges),'+','Color', Colors(iPlot,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
    text(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(2,Nges), ...
        PlanetPos(iPlot).xyz(3,Nges),Planets(iPlot).Name+'@T2', ...
        'VerticalAlignment','bottom','FontSize',8,'Color',Colors(iPlot,:))
end

x=linspace(PlanetPos(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe),10);
y=linspace(PlanetPos(3).xyz(2,ErdNaehe),Komet.xyz(2,ErdNaehe),10);
z=linspace(PlanetPos(3).xyz(3,ErdNaehe),Komet.xyz(3,ErdNaehe),10);
plot3(x,y,z,'Color',Colors(3,:),'LineWidth',2);
plot3(PlanetPos(3).xyz(1,ErdNaehe),PlanetPos(3).xyz(2,ErdNaehe), ...
    PlanetPos(3).xyz(3,ErdNaehe),'o','Color', Colors(AK,:), ...
    'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
text(PlanetPos(3).xyz(1,ErdNaehe),PlanetPos(3).xyz(2,ErdNaehe), ...
    PlanetPos(3).xyz(3,ErdNaehe),'Erdnaehe','VerticalAlignment', ...
    'bottom','FontSize',8,'Color',Colors(3,:));

% Einstellungen der Achsen 
xlim([-Ax_L Ax_L]);
ylim([-Ax_L Ax_L]);
zlim([-Ax_L Ax_L]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));

title(strcat(BaPaK.Name(AK),' - uebersicht'));
shading flat
alpha 0.1
grid on;
grid minor;
axis square;
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');

%%_________________________________________________________________________
% 
% Figure 3: Aequatoriale Koordinaten
figure('Name',header1);

RA(:)=wrapTo180(rad2deg(Komet.equ(2,:)));       % Rektaszension 
DEC(:)=rad2deg(Komet.equ(3,:));                 % Deklination 
plot(RA,DEC,'Color', Colors(AK,:));
hold on;
plot(RA,DEC,':+','MarkerIndices',1:NrStepsMarker: ...
    length(Komet.xyz(1,:)),'LineWidth',1,'Color',Colors(AK,:));

for k=1:MaxStepsLabel     
    dt4 = datetime(Komet.Time((k-1)*NrStepsLabel+1),'ConvertFrom', ...
        'juliandate');
    mylabels(k,:)= string(dt4,'dd-MM-yy');
end

for k=1:MaxStepsLabel
    Xlab(k)= RA((k-1)*NrStepsLabel+1);
    Ylab(k)=DEC((k-1)*NrStepsLabel+1);
end

LabelPoints(Xlab,Ylab,mylabels,'E',0.01,1,'FontSize',8,'Color', ...
    Colors(AK,:));
title(strcat(BaPaK.Name(AK),' - aequatoriale Koordinaten'));
hold on;
grid on;
grid minor;
xlabel('Rektaszension in °')
ylabel('Deklination in °');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
