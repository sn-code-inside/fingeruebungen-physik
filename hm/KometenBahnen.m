% -------------------------------------------------------------------------
% KometenBahnen.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahnen verschiedener Kometen in heliozentrischen
% Koordinaten unter Benutzung der Gaussschen Vektoren.
% (Aequinoktium und Ekliptik J2000) 
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
% Laden der Bahnparameter der Planeten und deren Ableitung
[BaPa,BaPadot]=OrbitParameter(T1, Aequi);
% Laden der Bahnparameter der Kometen 
BaPaK=KometParameter('Data/Kometen.csv');

AK=4; %Zahl der dargestellten Kometen

% -------------------------------------------------------------------------
% Beginn Rechnung

% Berechnung der Kometenbahnen nach Ellipsengleichung
Komets = KometOrbits(BaPaK);
% Berechnung der Planetenbahnen nach Ellipsengleichung
Planets = PlanetOrbits(BaPa);

% -------------------------------------------------------------------------
%%
% Graphische Ausgabe

% Figure 1: Kometenbahnen 
header1='Kometen: Bahnen in der Ekliptik';
figure('Name',header1);
subplot(1,2,1);         % linkes Schaubild in Figure 1 
% Iteration durch die Kometen 
for iPlot = 1:AK
    p(iPlot)=plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:), ...
        'Color',Colors(iPlot,:));
    p(iPlot).LineWidth=1;
    hold on
    axis equal
end
% Achse Sonne-Fruehlingspunkt 
SonneFP; 
% Iteration durch die Planeten 
for iPlot = 3:8
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
        Colors(iPlot,:));
    hold on;
    axis equal;
end
% Beschriftung der aeusseren Planetenbahnen 
for iPlot = 5:8
    h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1),'+', ...
        'C','FontSize',10,'Color',Colors(iPlot,:));
    h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
        Planets(iPlot).Name,'E','FontSize',8,'Color',Colors(iPlot,:));
    hold on;
    axis equal;
end

% Einstellung der Achsen 
ylim([-40 40])
xlim([-40 40]);
grid on;
grid minor;
% Beschriftungen 
header2 = strcat('Bahnen von Kometen');
title(header2);
xlabel('x in AE')
ylabel('y in AE')
% Iteration durch die Kometen fuer die Legende 
for iPlot=1:AK
    lgdtxt(iPlot) = Komets(iPlot).Name;
end
legend(lgdtxt,'location','northeast'),
legend boxoff;

subplot(1,2,2);         % rechtes Schaubild in Figure 1 
% Iteration durch die Kometen 
for iPlot = 1:AK
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(3,:),'Color', ...
        Colors(iPlot,:));
    hold on
   axis equal
end
% Achse Sonne-Fruehlingspunkt 
SonneFP;
% Einstellungen der Achsen 
xlim([-40 40])
ylim([-40 40]);
grid on;
grid minor;
% Beschriftungen 
title(header2);
xlabel('x in AE')
ylabel('z in AE');
legend(lgdtxt,'location','northeast'),
legend boxoff;

% Figure 2: Graphische Darstellung Halley 1986 und 1759

header1='Kometen im Vergleich';
figure('Name',header1);
subplot(1,2,1);         % linkes Schaubild in Figure 2 
% Iteration durch beide Halley-Bahnen 
for iPlot = 7:8
    p(iPlot)=plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:), ...
        'Color',Colors(iPlot,:));
    p(iPlot).LineWidth=1;
    hold on
    axis equal
end
% Achse Sonne-Fruehlingspunkt 
SonneFP;
% Iteration durch die Planeten 
for iPlot= 3:9
      plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:), ...
          Planets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
      hold on;
      axis equal;
end
% Beschriftung der Planetenbahnen 
for iPlot = 3:7
      h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
          '+','C','FontSize',10,'Color',Colors(iPlot,:));
      h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
          Planets(iPlot).Name,'E','FontSize',8,'Color',Colors(iPlot,:));
end
h=LabelPoints(-25,17,Planets(8).Name,'E','FontSize',8,'Color',Colors(8,:));
h=LabelPoints(-25,27,Planets(9).Name,'E','FontSize',8,'Color',Colors(9,:));

% Einstellungen der Achsen 
xlim([-30 5]);
ylim([-5 30])
grid on;
grid minor;
% Beschriftungen 
header2 = strcat('Halley Bahnen im Vergleich 1986 und 1758');
title(header2);
xlabel('x in AE')
ylabel('y in AE')
legend(Komets(7).Name,Komets(8).Name,'Sonne','Richtung aequinoktium', ...
    'location','northeast');
legend boxoff;

subplot(1,2,2);         % rechtes Schaubild in Figure 2 
% Iteration durch beide Halley-Bahnen 
for iPlot = 7:8
    plot(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(3,:),'Color', ...
        Colors(iPlot,:),'LineWidth',1);    
    hold on;
    axis equal;
end
% Achse Sonne-Fruehlingspunkt 
SonneFP;
% Einstellungen der Achsen 
xlim([-10 2])
ylim([-10 2]);
grid on;
grid minor;
% Beschriftungen 
title(header2);
xlabel('x in AE')
ylabel('z in AE');
legend(Komets(7).Name,Komets(8).Name,'Sonne','Richtung aequinoktium', ...
    'location','northeast');
legend boxoff;

%__________________________________________________________________________
%%
% Figure 3: 3-Dimensionale Darstellung
%
header2=' Bahnen von Kometen';
figure('Name',header2);
for resolution=2:2
    % Iteration durch die Kometen 
    for iPlot = 1:AK
        p=plot3(Komets(iPlot).xyz(1,:),Komets(iPlot).xyz(2,:), ...
            Komets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
        p.LineWidth=2;
        hold on
        axis equal
    end
    % Achse Sonne-Fruehlingspunkt 
    SonneFP;
    % Iteration durch die Planeten fuer die Planetenbahnen 
    for iPlot= 5-(resolution-1)*4:8-(resolution-1)*4
      plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:), ...
          Planets(iPlot).xyz(3,:),'Color', Colors(iPlot,:));
      hold on;
      axis equal;
    end
    % Iteration durch die Planeten fuer die Beschriftung 
    for iPlot = 5-(resolution-1)*4:8-(resolution-1)*4
        h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
            '+','C','FontSize',10,'Color',Colors(iPlot,:));
        h=LabelPoints(Planets(iPlot).xyz(1,1), Planets(iPlot).xyz(2,1), ...
            Planets(iPlot).Name,'E','FontSize',8,'Color',Colors(iPlot,:));
    end
    if resolution==1
        plot3(Planets(9).xyz(1,:),Planets(9).xyz(2,:), ...
            Planets(9).xyz(3,:),'Color', Colors(9,:));
        h=LabelPoints(Planets(9).xyz(1,91), Planets(9).xyz(2,91), ...
            Planets(9).Name,'E','FontSize',8,'Color',Colors(9,:));
    end
    % Einstellungen der Achsen 
    xl = xlim;
    yl = ylim;
    [X,Y] = meshgrid(xl,yl);
    surf(X,Y,zeros(size(X)))
    shading flat
    alpha 0.1
    grid on;
    grid minor;
    axis equal;
    if resolution ==1 
        xlim([-100 50]);
        ylim([-100 50]);
        zlim([-100 50]);
        title(strcat(header1,' - uebersicht'));
    else
        xlim([-3 2]);
        ylim([-3 2]);
        zlim([-1 1]);
        title(strcat(header1,' - Perihelumgebung'));
    end
     xlabel('x in AE')
    ylabel('y in AE');
    zlabel('z in AE');
end
txt = "";
legend(lgdtxt,'location','southeast'),
legend boxoff;
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
