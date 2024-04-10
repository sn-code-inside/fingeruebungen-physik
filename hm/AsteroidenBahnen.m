% -------------------------------------------------------------------------
% AsteroidenBahnen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahnen verschiedener Asteroiden in 
% heliozentrischen Koordinaten unter Benutzung der Gaussschen Vektoren.
% (Aequinoktium und Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung Zeit und Aequinoktium
T1 = juliandate(datetime(2000,1,1,12,0,0));
Aequi = 'J2000';

% Einlesen der Bahnparameter der Planeten und deren Ableitung 
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);

% Einlesen der Bahnparameter der Kometen 
BaPaK=KometParameter('Kometen.csv');

% Zahl der dargestellten Asteroiden 9 ... AK
nAst = 3; 
if (nAst > 0) && (nAst <=4)
    AK   = nAst + 8;
else
    AK = 11; 
end

% Darstellungsbereich in AE eingeben
AE = 4;

%-------------------------------------------------------------------------
%%
% Beginn Rechnung
% Berechnung der Asteroiden- und Kometenbahnen (Ellipsengleichung)

Komets = KometOrbits(BaPaK);

%  Berechnung der Planetenbahnen (Ellipsengleichung)
Planets=PlanetOrbits(BaPa);

%--------------------------------------------------------------------------
%%
% Graphische Ausgabe

% Figure 1: 2-Dimensionale Darstellung

header1='Asteroiden: Bahnen in der Ekliptik';
figure('Name',header1);

% x-y-Ebene
subplot(1,2,1);     % linkes Schaubild in Figure 1 

% Iteration durch die Asteroiden 
for iPlot = 9:AK
    p(iPlot)=plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:), ...
        'Color', Colors(iPlot,:));
    p(iPlot).LineWidth=1;
    hold on
    axis equal
end

% Achse Sonne-Fruehlingspunkt 
SonneFP;

% Iteration durch die Planeten fuer die Planetenbahnen 
for iPlot = 2:4
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ... 
        Colors(iPlot,:));
    hold on;
    axis equal;
end
% Iteration durch die Planeten fuer die Beschriftung 
for iPlot = 2:4
    h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
        '+','C','FontSize',10,'Color',Colors(iPlot,:));
    h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
        Planets(iPlot).Name,'E','FontSize',8,'Color',Colors(iPlot,:));
    hold on;
    axis equal;
end

% Einstellung der Achsen 
ylim([-AE AE])
xlim([-AE AE]);
grid on;
grid minor;
% Beschriftungen 
header2 = strcat('Bahnen von Asteroiden');
title(header2);
xlabel('x in AE')
ylabel('y in AE')
legend(Komets(9).Name,Komets(10).Name,Komets(11).Name,'Sonne', ... 
    'Richtung FP','location','northeast');
legend boxoff;

% x-z-Ebene
subplot(1,2,2);         % rechtes Schaubild in Figure 1
% Iteration durch die Asteroiden 
for iPlot = 9:AK
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(3,:),'Color', ... 
        Colors(iPlot,:),'LineWidth',1);
    hold on
   axis equal
end

% Achse Sonne-Fruehlingspunkt 
SonneFP;

% Einstellung der Achsen
xlim([-AE AE])
ylim([-AE AE]);
grid on;
grid minor;
% Beschriftungen 
title(header2);
xlabel('x in AE')
ylabel('z in AE');
legend(Komets(9).Name,Komets(10).Name,Komets(11).Name,'Sonne', ... 
    'Ebene der Ekliptik','location','northeast');
legend boxoff;

% Figure 2: 3-Dimensionale Darstellung
figure('Name',header1);

% Iteration durch die Asteroiden 
for iPlot= 9:AK
  p=plot3(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:), ... 
      Komets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
  p.LineWidth=2;
  hold on;
  axis square;
end
hold on
axis equal

% Achse Sonne-Fruehlingspunkt
SonneFP,

% Iteration durch die Planeten fuer die Planetenbahnen
for iPlot= 3:4
  plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:), ...
      Planets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
  hold on;
  axis square;
end
% Iteration durch die Planeten fuer die Beschriftung 
for iPlot = 3:4
    text(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
        Planets(iPlot).xyz(3,1),'<','FontSize',8,'Color',Colors(iPlot,:))
    text(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
        Planets(iPlot).xyz(3,1),Planets(iPlot).Name, ...
        'VerticalAlignment','bottom','FontSize',8,'Color',Colors(iPlot,:));
end
hold on;
axis equal;

% Einstellungen der Achsen 
xlim([-AE AE]);
ylim([-AE AE]);
% zlim([-AE AE]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)));

% Beschriftungen 
title(strcat('Asteroiden - Uebersicht'));
legend(Komets(9).Name,Komets(10).Name,Komets(11).Name,'location', ...
    'northeast');
legend boxoff;
shading flat
alpha 0.1
grid on;
grid minor;
axis square;
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in AE');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
