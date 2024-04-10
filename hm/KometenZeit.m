% -------------------------------------------------------------------------
% KometenZeit.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahn eines Kometen als Funktion der Zeit 
% unter fallweiser Benutzung einer parabolischen Naeherung.
% (Aequinoktium und Ekliptik Datum bzw. J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
% Einlesen der Bahnparameter der Kometen und Asteroiden 
BaPaK=KometParameter('Data/Kometen.csv'); 
GM      = 2.95479E-04;     % µG [AE^3/d^2]
c_l     = 173.14;          % Lichtgeschwindigkeit [AE/d]

% Initialisierung & Input

% Input
% Auswahl des Komet (Nr 1-7 eingeben)
HK=2;
% Zeitraum und Darstellung Perihelnaehe Neowise
% fuer andere Kometen ggf. aendern (siehe unten) 
dt1 = datetime('2020-03-01 00:00:00');
dt2 = datetime('2020-09-01 00:00:00');
Aequi = 'Datum';

% Zeitraum und Darstellung ueber > 10 Jahre Neowise
% HK=2;
% dt1 = datetime('2005-02-01 00:00:00');
% dt2 = datetime('2045-02-01 00:00:00');
% Aequi = 'J2000';

% Zeitraum und Darstellung Perihelnaehe Halley
% HK=1;
% dt1 = datetime('1985-11-01 00:00:00');
% dt2 = datetime('1986-05-01 00:00:00');
% Aequi = 'Datum';

% oder
% Zeitraum und Darstellung ueber > 10 Jahre Halley
% HK=1;
% dt1 = datetime('1985-02-01 00:00:00');
% dt2 = datetime('2063-02-01 00:00:00');
% Aequi = 'J2000';

T1 = juliandate(dt1);
T2 = juliandate(dt2);
% Einlesen der Bahnparameter der Planeten und deren Ableitungen 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi); 

if T2-T1 <= 700
    % Berechnung fuer kurze Zeitraeume 
    long=0;
    Nges=10000; %Zahl der Stuetzstellen im Intervall T1 - T2
    T=linspace(T1,T2,Nges); %Berechnungsintervall
    StepsPerDay= Nges/(T2-T1); %Zahl der Stuetzstellen pro Tag
    DatumStr1 =  string(dt1,'dd.MM.yy');
    DatumStr2 =  string(dt2,'dd.MM.yy');
    %Input
    StepsMarker= 5; %Zahl der markierten Tage 
    %(Richtwert 5 bei T2-T1=1 Jahr usw.)
    NrStepsMarker=round(StepsMarker*StepsPerDay); %Zahl der markierten Tage
    %Input
    StepsLabel=2; %Zahl der angezeigten Tage
    %((Richtwert 2 T2-T1=1 Jahr usw.)
    NrStepsLabel=StepsLabel*NrStepsMarker; %Zahl der angezeigten Tage
    MaxStepsLabel=floor(Nges/NrStepsLabel);
    AE=2.5;
    ia=3;
    ie=4;
else
    % Berechnung fuer lange Zeitraeume 
    long=1;
    DatumStr1 =  string(dt1,'yyyy');
    DatumStr2 =  string(dt2,'yyyy');
    %Input
    StepsperDay = 1;
    StepsperYear = 365.25;
    Nges = round(T2-T1); %Zahl der Stuetzstellen im Intervall T1 - T2
    StepsMarker= 2; %Zahl der markierten Jahre 
    %(Richtwert 2 bei T2-T1=100 Jahre usw.)
    NrStepsMarker=round(StepsMarker*StepsperYear); %Anzahl markierte Tage
    StepsLabel=5; %Zahl der beschrifteten Jahre
    %((Richtwert 5 T2-T1=100 Jahre usw.)
    NrStepsLabel=StepsLabel*NrStepsMarker; %Zahl der angezeigten Jahre
    MaxStepsLabel=floor(Nges/NrStepsLabel)+1;
    T=linspace(T1,T1+Nges,Nges); %Berechnungsintervall
    AE=40;
    ia=5;
    ie=8;
end
T_Per =BaPaK.T0(HK);       % Perihelzeit des Kometen


%-------------------------------------------------------------------------
%%
%Beginn Rechnung

%Berechnung der Erdposition zum Perihel des Kometen
ErdePeri(3)=PlanetPQR(T_Per, BaPa, BaPadot, 3);

% Berechnung der heliozentrischen Koordinaten der Erde und 
% der anderen Planeten von T1 bis T2 
for k=1:8
    PlanetPos(k)=PlanetPQR(T, BaPa, BaPadot,k);
end

% Berechnung der heliozentrischen Koordinaten 
% des Kometen nach Loesung mit Stumppfschen Funktionen
% von T1 bis T2 

Komet  = KometPQR(GM, T, BaPaK,HK);

%PQR Multiplikation Perihel
yin=[BaPaK.qP(HK);0;0];
OmegaPn=deg2rad(BaPaK.OmegaP(HK));
omegaPn=deg2rad(BaPaK.omegaP(HK));
iPn=deg2rad(BaPaK.iP(HK));
Perihel(HK).xyz=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);

% Berechnung der Planetenbahnen Ellipsengleichung
Planets=PlanetOrbits(BaPa);

%Berechnung Abstand Erde-Komet
rek=Komet.xyz-PlanetPos(3).xyz;
rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));
[minrekbetr,ErdNaehe] = min(rekbetr);

%Beruecksichtigung Lichtlaufzeit
rek_korr= rek-rekbetr.*Komet.v/c_l;

%aequatoriale Koordinaten
eps0 = deg2rad(EpsErde(T1));
% Koordinatentrafo in Equ-System
Komet.equ = CalcAnglesfromXYZ(mtimes(R_x(-eps0),rek_korr)); 

%-------------------------------------------------------------------------
%% 
% Graphische Ausgabe

% Figure 1: 2-D-Ausgabe Kometenbahnen 
header1='Kometbahnen';
figure('Name',header1);

% x-y-Ebene
subplot(1,2,1);         % linkes Schaubild in Figure 1 
% Plot Komet
plot(Komet.xyz(1,:),Komet.xyz(2,:),':+','MarkerIndices', ...
    1:NrStepsMarker:length(Komet.xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(HK,:));
hold on;
% Einstellungen der Achsen 
xlim([-AE AE]);
ylim([-AE AE]);
grid on;
grid minor;
% Beschriftungen der Zeitpunkte 
for k=1:MaxStepsLabel
    dt4 = datetime(Komet.Time((k-1)*NrStepsLabel+1),'ConvertFrom', ...
        'juliandate');
    if long == 0 
        mylabels(k,:)= string(dt4,'dd-MM-yy');
    else
        mylabels(k,:)= string(dt4,'yyyy');
    end
end
for k=1:MaxStepsLabel
    Xlab(k)=Komet.xyz(1,(k-1)*NrStepsLabel+1);
    Ylab(k)=Komet.xyz(2,(k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,Ylab,mylabels,'E',0.01,0,'FontSize',8,'Color', ...
    Colors(HK,:));
 
% Sonne, Perihel , Richtung Fruehlingspunkt
SonneFP;
x=linspace(0,Perihel(HK).xyz(1),10);
y=linspace(0,Perihel(HK).xyz(2),10);
p(10)=plot(x,y,'Color',Colors(10,:));
p(10).LineWidth=2;

% Iteration durch die Planeten fuer Planetenbahnen
for iPlot = 2:8
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
        Colors(iPlot,:));
end

if long == 0 
    % fuer kurze  Zeitraeume
    % Position Erdnaehe
    x=linspace(PlanetPos(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe),10);
    y=linspace(PlanetPos(3).xyz(2,ErdNaehe),Komet.xyz(2,ErdNaehe),10);
    p(10)=plot(x,y,'Color',Colors(3,:));
    p(10).LineWidth=2;
    plot(PlanetPos(3).xyz(1,ErdNaehe),PlanetPos(3).xyz(2,ErdNaehe),'o', ...
        'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(3,:));
    h=LabelPoints(PlanetPos(3).xyz(1,ErdNaehe), ...
        PlanetPos(3).xyz(2,ErdNaehe),'Erdnaehe','E',0.02,'FontSize',8, ...
        'Color',Colors(3,:));
    % Position der Planeten Mars und Erde am Anfang und Ende des Intervalls
    plot(ErdePeri(3).xyz(1,1),ErdePeri(3).xyz(2,1),'+','Color', ...
        Colors(3,:),'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
    h=LabelPoints(ErdePeri(3).xyz(1,1), ErdePeri(3).xyz(2,1), ...
        'Erde@Komet-Perihel','S',0.02,'FontSize',8,'Color',Colors(3,:));
    % Planetenbahnen Mars und Erde 
    for iPlot = 3:4
        plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1),'+', ...
            'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
            Colors(iPlot,:));
        h=LabelPoints(PlanetPos(iPlot).xyz(1,1), ...
            PlanetPos(iPlot).xyz(2,1),Planets(iPlot).Name+'@T1','E', ...
            0.02,'FontSize',8,'Color',Colors(iPlot,:));
        plot(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(2,Nges), ...
            '+','Color', Colors(iPlot,:),'MarkerSize',5, ...
            'MarkerFaceColor',Colors(iPlot,:));
        h=LabelPoints(PlanetPos(iPlot).xyz(1,Nges), ...
            PlanetPos(iPlot).xyz(2,Nges),Planets(iPlot).Name+'@T2','E', ...
            0.02,'FontSize',8,'Color',Colors(iPlot,:));
    end
end
% Beschriftungen 
header2 = strcat(BaPaK.Name(HK));
title(header2);
axis square;
xlabel('x in AE');
ylabel('y in AE');

%________________________________________________________________________

% x-z-Ebene
subplot(1,2,2);         % rechtes Schaubild Figure 1 
p(HK)=plot(Komet.xyz(1,:),Komet.xyz(3,:),'Color', Colors(HK,:));
% Plot Komet
hold on;
plot(Komet.xyz(1,:),Komet.xyz(3,:),'-+','MarkerIndices', ...
    1:NrStepsMarker:length(Komet.xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(HK,:));
axis square;
grid on;
grid minor;

% Beschriftungen der Zeitpunkte 
for k=1:MaxStepsLabel     
    dt4 = datetime(Komet.Time((k-1)*NrStepsLabel+1),'ConvertFrom', ...
        'juliandate');
    if long == 0 
        mylabels(k,:)= string(dt4,'dd-MM-yy');
    else
        mylabels(k,:)= string(dt4,'yyyy');
    end
end
for k=1:MaxStepsLabel
    Xlab(k)=Komet.xyz(1,(k-1)*NrStepsLabel+1);
    Ylab(k)=Komet.xyz(3,(k-1)*NrStepsLabel+1);
end
LabelPoints(Xlab,Ylab,mylabels,'E',0.01,0,'FontSize',8,'Color', ...
    Colors(HK,:));

% Sonne, Perihel , Richtung Fruehlingspunkt
SonneFP;
x=linspace(0,Perihel(HK).xyz(1),10);
y=linspace(0,Perihel(HK).xyz(3),10);
p(10)=plot(x,y,'Color',Colors(10,:));
p(10).LineWidth=2;

% Planetenbahnen
if long == 0 
    iende = 4; % Erde und Mars fuer kurze Zeitraeume
else
    iende = 8; % weitere Planeten fuer lange Zeitraeume 
end
for iPlot = 3:iende
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(3,:),'Color', ...
        Colors(iPlot,:));
end

if long == 0
    % Berechnung fuer kurze Zeitraeume 
    % Position Erdnaehe
    x=linspace(PlanetPos(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe),10);
    z=linspace(PlanetPos(3).xyz(3,ErdNaehe),Komet.xyz(3,ErdNaehe),10);
    p(10)=plot(x,z,'Color',Colors(3,:));
    p(10).LineWidth=2;
    plot(PlanetPos(3).xyz(1,ErdNaehe),PlanetPos(3).xyz(3,ErdNaehe),'o', ...
        'Color',Colors(HK,:),'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
    h=LabelPoints(PlanetPos(3).xyz(1,ErdNaehe), ...
        PlanetPos(3).xyz(3,ErdNaehe),'Erdnaehe','E',0.02,'FontSize',8, ...
        'Color',Colors(3,:));
    % Position Planeten Anfang und Ende Intervall, Erde am Perihel
    plot(ErdePeri(3).xyz(1,1),ErdePeri(3).xyz(3,1),'+','Color', ...
        Colors(3,:),'MarkerSize',5,'MarkerFaceColor',Colors(3,:));
    h=LabelPoints(ErdePeri(3).xyz(1,1), ErdePeri(3).xyz(3,1), ...
        'Erde@Komet-Perihel','S',0.02,'FontSize',8,'Color',Colors(3,:));
    % Beschriftung der Planetenbahnen 
    for iPlot = 3:iende 
        plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(3,1),'+', ...
            'Color', Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
            Colors(iPlot,:));
        h=LabelPoints(PlanetPos(iPlot).xyz(1,1), ...
            PlanetPos(iPlot).xyz(3,1),Planets(iPlot).Name+'@T1','E', ...
            0.02,'FontSize',8,'Color',Colors(iPlot,:));
        plot(PlanetPos(iPlot).xyz(1,Nges),PlanetPos(iPlot).xyz(3,Nges), ...
            '+','Color', Colors(iPlot,:),'MarkerSize',5, ...
            'MarkerFaceColor',Colors(iPlot,:));
        h=LabelPoints(PlanetPos(iPlot).xyz(1,Nges), ...
            PlanetPos(iPlot).xyz(3,Nges),Planets(iPlot).Name+'@T2','E', ...
            0.02,'FontSize',8,'Color',Colors(iPlot,:));
    end
end

% Beschriftung 
header2 = strcat(BaPaK.Name(HK));
title(header2);
xlim([-AE AE]);
% ylim([-AE AE]);
xlabel('x in AE');
ylabel('z in AE');

%_________________________________________________________________

% Figure 2: 3-Dimensionale Darstellung
figure('Name',header1);

for resolution=1:2 % 1 Übersicht % 2 Um Perihel herum
    subplot(1,2,resolution);
    p=plot3(Komet.xyz(1,:),Komet.xyz(2,:),Komet.xyz(3,:),'Color', ...
        Colors(HK,:));
    p.LineWidth=2;
    hold on
    axis equal
    % Achse Sonne-Fruehlingspunkt 
    SonneFP;
    for iPlot= 5-(resolution-1)*4:8-(resolution-1)*4
      plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:), ...
          Planets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
      hold on;
      axis square;
    end
    for iPlot = 5-(resolution-1)*4:8-(resolution-1)*4
        plot3(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1), ...
            PlanetPos(iPlot).xyz(3,1),'o','Color', Colors(iPlot,:), ...
            'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
        h=LabelPoints(PlanetPos(iPlot).xyz(1,1), ...
            PlanetPos(iPlot).xyz(2,1),Planets(iPlot).Name,'S',0.03, ...
            'FontSize',8,'Color',Colors(iPlot,:));
    end
   x=linspace(PlanetPos(3).xyz(1,ErdNaehe),Komet.xyz(1,ErdNaehe),10);
   y=linspace(PlanetPos(3).xyz(2,ErdNaehe),Komet.xyz(2,ErdNaehe),10);
   z=linspace(PlanetPos(3).xyz(3,ErdNaehe),Komet.xyz(3,ErdNaehe),10);
   if resolution==2 
        plot3(x,y,z,'Color',Colors(3,:),'LineWidth',2);
        plot3(PlanetPos(3).xyz(1,ErdNaehe), ...
            PlanetPos(3).xyz(2,ErdNaehe),PlanetPos(3).xyz(3,ErdNaehe), ...
            'o','Color', Colors(iPlot,:),'MarkerSize',5, ...
            'MarkerFaceColor',Colors(3,:));
        h=LabelPoints(PlanetPos(3).xyz(1,ErdNaehe), ...
            PlanetPos(3).xyz(2,ErdNaehe),Planets(3).Name,'E',0.01, ...
            'FontSize',8,'Color',Colors(3,:));    
        x=linspace(0,Perihel(HK).xyz(1),10);
        y=linspace(0,Perihel(HK).xyz(2),10);
        z=linspace(0,Perihel(HK).xyz(3),10);
        p(10)=plot3(x,y,z,'Color',Colors(10,:));
        p(10).LineWidth=2;
   end
   % Einstellungen der Achsen 
   xlim([-40 40]);
   ylim([-40 40]);
   xl = xlim;
   yl = ylim;
   [X,Y] = meshgrid(xl,yl);
   surf(X,Y,zeros(size(X)));
   zlim([-40 40]);
   shading flat
   alpha 0.1
   grid on;
   grid minor;
   axis square;
   if resolution ==1 
        AEax = 40;
        title(strcat(BaPaK.Name(HK),' - uebersicht'));
   else
        AEax = 2;
        title(strcat(BaPaK.Name(HK),' - Perihelumgebung'));
   end
   xlim([-AEax AEax]);
   ylim([-AEax AEax]);
   zlim([-AEax AEax]);
   xlabel('x in AE')
   ylabel('y in AE');
   zlabel('z in AE');
 end


%_________________________________________________________________________
% 
% Figure 3: Aequatoriale Koordinaten
 
if long == 0
    % Berechnung fuer kurze Zeitraeume 
    figure('Name',header1);
    RA(:)=wrapTo360(rad2deg(Komet.equ(2,:)));
    DEC(:)=rad2deg(Komet.equ(3,:));
    plot(RA,DEC,'Color', Colors(HK,:));
    hold on;
    plot(RA,DEC,':+','MarkerIndices',1:NrStepsMarker: ...
        length(Komet.xyz(1,:)),'LineWidth',1,'Color',Colors(HK,:));
    % Beschriftung der Zeitpunkte 
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
        Colors(HK,:));
    title(strcat(BaPaK.Name(HK),'-aequatoriale Koordinaten'));
    hold on;
    grid on;
    grid minor;
    xlabel('Rektaszension in °')
    ylabel('Deklination in °');
    set(gca,'FontSize',16);

end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
