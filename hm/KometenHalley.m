% -------------------------------------------------------------------------
% KometenHalley.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahnen des Kometen Halley von 1986 und 1758
% auf Basis der Keplergleichung.
% (Aequinoktium und die Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
T1 = juliandate(datetime(2000,1,1,12,0,0));
Aequi = 'J2000';
% Einlesen der Bahnparameter der Planeten und deren Ableitung 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);
% Einlesen der Bahnparameter der Kometen 
BaPaK=KometParameter('Kometen.csv');

GM      = 2.95479E-04;     % µG [AE^3/d^2]
c0      = 0.00578;         % Lichtgeschwindigkeit 
As      = 2.5;             % Limits der Achsen 

%Berechnung der Planetenbahnen nach Keplerloesung
Planets=PlanetOrbits(BaPa);

%-------------------------------------------------------------------------
%%
% Beginn Rechnung
% Auswahl des Kometen 
HK=6;
iPlot = 1;

while iPlot < 3    
    % Auswahl der Kometen Halley aus der Liste
    HK=6+iPlot;
    name(iPlot,:)=BaPaK.Name(HK);
    iPlot=iPlot+1;
    if HK==7
        dt1 = datetime('1985-07-01 00:00:00'); 
        T1 = juliandate(dt1); % Julianisches Datum
        dt2 = datetime('1986-06-30 00:00:00'); 
        T2 = juliandate(dt2); % Julianisches Datum
    else
        dt1 = datetime('1758-07-01 00:00:00'); 
        T1 = juliandate(dt1); % Julianisches Datum
        dt2 = datetime('1759-06-30 00:00:00'); 
        T2 = juliandate(dt2); % Julianisches Datum
    end   
    T0    = (T1-2451545)/36525;
    eps0  = deg2rad(EpsErde(T0)); % Schiefe der Ekliptik 
    Nges=10000;
    T=linspace(T1,T2,Nges)+1;
    % Einstellen der Schrittweiten 
    StepsperDay=round(Nges/(T2-T1));
    NrStepsMarker=5*StepsperDay;
    NrStepsLabel=NrStepsMarker*2;
    MaxStepsLabel=floor(Nges/NrStepsLabel);
    
    T_Per(HK)=BaPaK.T0(HK);
    % Berechnung der Planetenposition zum Perihel, zu T1 und zu T2
    % Iteration durch die Planeten 
    for k=1:8
        PlanetPos(k)=PlanetPQR(T_Per(HK), BaPa, BaPadot, k);
    end

    % Berechnung der heliozentrischen Koordinaten der Erde und 
    % des Kometen nach Keplerloesung
    % von T1 bis T2 
    ePn=BaPaK.eP(HK);
    if HK==7 
        Halley1986=KometPQR(GM, T, BaPaK,HK);         
        Komets=[Halley1986, Halley1986];
        Erde1=PlanetPQR(T, BaPa, BaPadot, 3);
        Erde=[Erde1, Erde1];
     else
        Halley1758=KometPQR(GM, T, BaPaK,HK);
        Komets=[Halley1986, Halley1758];
        Erde2=PlanetPQR(T, BaPa, BaPadot, 3);
        Erde=[Erde1, Erde2];
    end
end

for iPlot=1:2
    HK=6+iPlot;
    OmegaPn=deg2rad(BaPaK.OmegaP(HK));
    omegaPn=deg2rad(BaPaK.omegaP(HK));
    iPn=deg2rad(BaPaK.iP(HK));
    % PQR Multiplikation Perihel
    yin=[BaPaK.qP(HK);0;0];
    Perihel(iPlot).xyz=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
    % Berechnung Abstand Erde-Komet
    rek=Komets(iPlot).xyz-Erde(iPlot).xyz;
    rekbetr=sqrt(rek(1,:).*rek(1,:)+rek(2,:).*rek(2,:)+rek(3,:).*rek(3,:));
    [minrekbetr,ErdNaehe(iPlot)] = min(rekbetr);
    % Beruecksichtigung Lichtlaufzeit
    rek_korr= rek-c0*rekbetr.*Komets(iPlot).v;
    % Koordinatentransformation in aequatoriale Koordinaten
    Komets(iPlot).equ = CalcAnglesfromXYZ(mtimes(R_x(-eps0),rek_korr)); 
end
%--------------------------------------------------------------------------
% Graphische Ausgabe

header1='Kometen: Bahnen in der Ekliptik';
figure('Name',header1);
% x-y-Ebene 
subplot(1,2,1);         % linkes Schaubild 

for iPlot =1:2
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:),':+', ...
        'MarkerIndices',1:NrStepsMarker:length(Komets(iPlot).xyz(1,:)), ...
        'LineWidth',1,'Color',Colors(iPlot,:));
    hold on;
end
plot(Erde(1).xyz(1,:),Erde(1).xyz(2,:),':+','MarkerIndices', ...
    1:NrStepsMarker*3:length(Erde(1).xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(3,:));
plot(Erde(2).xyz(1,:),Erde(2).xyz(2,:),':d','MarkerIndices', ...
    1:NrStepsMarker*3:length(Erde(2).xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(3,:));

% Achse Sonne-Fruehlingspunkt  
SonneFP;

for iPlot = 3:5
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
        Colors(iPlot,:));
end
 
for iPlot=1:2
    x=linspace(0,Perihel(iPlot).xyz(1),10);
    y=linspace(0,Perihel(iPlot).xyz(2),10);
    plot(x,y,'Color',Colors(10,:),'LineWidth',2);
    x=linspace(Erde(iPlot).xyz(1,ErdNaehe(iPlot)), ...
        Komets(iPlot).xyz(1,ErdNaehe(iPlot)),10);
    y=linspace(Erde(iPlot).xyz(2,ErdNaehe(iPlot)), ...
        Komets(iPlot).xyz(2,ErdNaehe(iPlot)),10);
    plot(x,y,'Color',Colors(HK,:),'LineWidth',2);
    plot(Erde(iPlot).xyz(1,ErdNaehe(iPlot)), ...
        Erde(iPlot).xyz(2,ErdNaehe(iPlot)),'o','Color',Colors(iPlot,:),...
        'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
    h=LabelPoints(Komets(iPlot).xyz(1,ErdNaehe(iPlot)), ...
        Komets(iPlot).xyz(2,ErdNaehe(iPlot)),strcat('Erdnaehe - ', ...
        name(iPlot)),'S',0.01,'FontSize',8,'Color',Colors(iPlot,:));
end

% Einstellungen der Achsen 
xlim([-As As]);
ylim([-As As]);
axis square;

ia      = 3;               
ie      = 4;

for iPlot = ia:ie
    plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(2,1),'o', ...
        'Color',Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,1), PlanetPos(iPlot).xyz(2,1),...
        Planets(iPlot).Name,'S',0.02,'FontSize',8,'Color',Colors(iPlot,:));
end

grid on;
grid minor;

% Beschriftungen 
header2 = strcat(BaPaK.Name(7),'-',BaPaK.Name(8),' im Vergleich ');
title(header2);
legend(string(BaPaK.Name(7)),string(BaPaK.Name(8)),'Erde 1986','Erde 1758');
legend box off;
xlabel('x in AE');
ylabel('y in AE');
set(gca,'FontSize',16);

%________________________________________________________________________
% x-z-Ebene
subplot(1,2,2);         % rechtes Schaubild 
for iPlot =1:2
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(3,:),':+', ...
        'MarkerIndices',1:NrStepsMarker:length(Komets(iPlot).xyz(1,:)), ...
        'LineWidth',1,'Color',Colors(iPlot,:));
    hold on;
end
plot(Erde(1).xyz(1,:),Erde(1).xyz(3,:),':+','MarkerIndices', ...
    1:NrStepsMarker*3:length(Erde(1).xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(3,:));
plot(Erde(2).xyz(1,:),Erde(2).xyz(3,:),':d','MarkerIndices', ...
    1:NrStepsMarker*3:length(Erde(2).xyz(1,:)),'LineWidth',1,'Color', ...
    Colors(3,:));

% Achse Sonne-Fruehlingspunkt  
SonneFP;

for iPlot=1:2
    x=linspace(0,Perihel(iPlot).xyz(1),10);
    y=linspace(0,Perihel(iPlot).xyz(3),10);
    plot(x,y,'Color',Colors(10,:),'LineWidth',2);
    x=linspace(Erde(iPlot).xyz(1,ErdNaehe(iPlot)), ...
        Komets(iPlot).xyz(1,ErdNaehe(iPlot)),10);
    y=linspace(Erde(iPlot).xyz(3,ErdNaehe(iPlot)), ...
        Komets(iPlot).xyz(3,ErdNaehe(iPlot)),10);
    p(10)=plot(x,y,'Color',Colors(HK,:));
    p(10).LineWidth=2;
    plot(Erde(iPlot).xyz(1,ErdNaehe(iPlot)), ...
        Erde(iPlot).xyz(3,ErdNaehe(iPlot)),'o','Color',Colors(iPlot,:), ...
        'MarkerSize',5,'MarkerFaceColor',Colors(iPlot,:));
    h=LabelPoints(Komets(iPlot).xyz(1,ErdNaehe(iPlot)), ...
        Komets(iPlot).xyz(3,ErdNaehe(iPlot)),strcat('Erdnaehe - ', ...
        name(iPlot)),'E',0.01,'FontSize',8,'Color',Colors(iPlot,:));
end

for iPlot = 3:5
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(3,:),'Color', ...
        Colors(iPlot,:));
 end
for iPlot = 3:5
    plot(PlanetPos(iPlot).xyz(1,1),PlanetPos(iPlot).xyz(3,1),'o', ...
        'Color',Colors(iPlot,:),'MarkerSize',5,'MarkerFaceColor', ...
        Colors(iPlot,:));
    h=LabelPoints(PlanetPos(iPlot).xyz(1,1), PlanetPos(iPlot).xyz(3,1), ...
        Planets(iPlot).Name,'S',0.02,'FontSize',8,'Color',Colors(iPlot,:));
end

% Einstellungen der Achsen 
xlim([-As As]);
ylim([-1 0.5]);
axis square;
grid on;
grid minor;
% Legende 
legend(string(BaPaK.Name(7)),string(BaPaK.Name(8)),'Erde 1986','Erde 1758');
legend box off;

header2 = strjoin([BaPaK.Name(7),' - ',BaPaK.Name(8),' im Vergleich ']);
title(header2);
xlabel('x in AE');
ylabel('z in AE');
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
