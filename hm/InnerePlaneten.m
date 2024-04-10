% -------------------------------------------------------------------------
% InnerePlaneten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bahnen der inneren Planeten im Jahr 2004
% auf Basis der Keplergleichung.
% (Aequinoktium und Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
dt1 = datetime('2004-01-01 06:00:00');
T1 = juliandate(dt1); % Julianisches Datum %Startdatum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

dt2 = datetime('2004-06-08 06:00:00');
T2 = juliandate(dt2); % Julianisches Datum  %Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

dt3 = datetime('2005-01-01 06:00:00');
T3 = juliandate(dt3); % Julianisches Datum %Endatum
MJuDa3 = juliandate(dt3,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr3 =  string(dt3,'dd.MM.yyyy');
Aequi = 'Datum';
% Einlesen der Bahnparameter und deren Ableitung
[BaPa,BaPadot]=OrbitParameter(T2,Aequi);

%-------------------------------------------------------------------------
%Beginn Rechnung
maxP=1000;
T_vector=linspace(T1,T2,maxP); %Von Start bis  Bedeckung
T_vector2=linspace(T2,T3,maxP); %Von Bedeckung bis Ende
% Berechnung nach Keplerloesung

for k=1:4
    Planets(k)=PlanetPQR(T_vector, BaPa, BaPadot, k);
    Planets2(k)=PlanetPQR(T_vector2, BaPa, BaPadot, k);
end
%------------------------------------------------------------------------------
% Graphische Ausgabe

header1='Innere Planeten: Projektion Bahnen auf die Ekliptik';
figure('Name',header1);
%Planetenbahnen
for iPlot = 1:size(Planets,2)
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
         Colors(iPlot,:));
    hold on
    axis equal
end
%Sonne Fruehlingspunkt
SonneFP;  
%Transitlinie
plot([0 1.2*Planets(3).xyz(1,maxP)], [0 1.2*Planets(3).xyz(2,maxP)],...
      ':','Color', Colors(3,:),'LineWidth',2);
%Positionen Anfangs, Endpunkte, Midpoint
for iPlot = 1:4
    p=plot(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1),...
           'Color', Colors(iPlot,:));
    p.Marker='+';
    p.LineWidth=1;
    p=plot(Planets(iPlot).xyz(1,maxP),Planets(iPlot).xyz(2,maxP),...
           'Color', Colors(iPlot,:));
    p.Marker='o';
    p.LineWidth=2;
    p=plot(Planets2(iPlot).xyz(1,maxP),Planets2(iPlot).xyz(2,maxP),...
          'Color', Colors(iPlot,:));
    p.Marker='d';
    p.LineWidth=1;
   
end
for iPlot = 1:size(Planets,2)
    plot(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:),'Color',...
         Colors(iPlot,:),'LineStyle','-', 'LineWidth',1);
end
%Linien Anfangs, Endpunkte, Midpoint
for iPlot = 1:size(Planets,2)
    plot([0 Planets(iPlot).xyz(1,1)], [0 Planets(iPlot).xyz(2,1)],':',...
          'Color', Colors(iPlot,:));
    plot([0 Planets(iPlot).xyz(1,maxP)], [0 Planets(iPlot).xyz(2,maxP)],...
          '--','Color', Colors(iPlot,:));
    plot([0 Planets2(iPlot).xyz(1,end)], [0 Planets2(iPlot).xyz(2,end)],...
          '-.','Color', Colors(iPlot,:));
 end
ylim([-1.5 2])
xlim([-2 2]);
grid on;
grid minor,
header2 = strjoin([DatumStr1,' - ',DatumStr3]);
title(header2);
xlabel('x in AE')
ylabel('y in AE');
legend(Planets.Name,'Sonne','Richtung Frühlingspunkt',...
       'Venus Transit','location','southeast');
legend boxoff;
set(gca,'FontSize',16);

% 3-Dimensionale Darstellung
header1=' Bahnen der Inneren Planeten';
figure('Name',header1);
for iPlot = 1:size(Planets,2)
    p=plot3(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),...
            Planets(iPlot).xyz(3,:)*10,'Color', Colors(iPlot,:));
    p.LineWidth=2;
    hold on
    axis equal
end
SonneFP;
plot3([0 1.5*Planets(3).xyz(1,maxP)], [0 1.5*Planets(3).xyz(2,maxP)],...
    [0 1.5*Planets(3).xyz(3,maxP)],':','Color', Colors(3,:),'LineWidth',2);
for iPlot = 1:size(Planets,2)
    p=plot3(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:),...
    Planets2(iPlot).xyz(3,:)*10,'Color', Colors(iPlot,:),'LineStyle','--');
    p.LineWidth=2;
    hold on
    axis equal
end
for iPlot = 1:4
    p=plot3(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1),...
        10*Planets(iPlot).xyz(3,1),'Color', Colors(iPlot,:));
    p.Marker='+';
    p.LineWidth=1;
    p=plot3(Planets(iPlot).xyz(1,maxP),Planets(iPlot).xyz(2,maxP),...
        10*Planets(iPlot).xyz(3,maxP),'Color', Colors(iPlot,:));
    p.Marker='o';
    p.LineWidth=2;
    p=plot3(Planets2(iPlot).xyz(1,maxP),Planets2(iPlot).xyz(2,maxP),...
        10*Planets2(iPlot).xyz(3,maxP),'Color', Colors(iPlot,:));
    p.Marker='d';
    p.LineWidth=1;
   
 end
for iPlot = 1:size(Planets,2)
    plot3(Planets2(iPlot).xyz(1,:),Planets2(iPlot).xyz(2,:),...
    10*Planets2(iPlot).xyz(3,:),'Color', Colors(iPlot,:),'LineStyle','--');
end
for iPlot = 1:size(Planets,2)
    plot3([0 Planets(iPlot).xyz(1,1)],[0 Planets(iPlot).xyz(2,1)],...
        [0 10*Planets(iPlot).xyz(3,1)],':','Color', Colors(iPlot,:));
    plot3([0 Planets(iPlot).xyz(1,maxP)],[0 Planets(iPlot).xyz(2,maxP)],...
        [0 10*Planets(iPlot).xyz(3,maxP)],'--','Color', Colors(iPlot,:));
    plot3([0 Planets2(iPlot).xyz(1,end)],[0 Planets2(iPlot).xyz(2,end)],...
        [0 10*Planets2(iPlot).xyz(3,end)],'-.','Color', Colors(iPlot,:));
 end
xlim([-2 2]);
ylim([-2 2]);
xl = xlim;
yl = ylim;
[X,Y] = meshgrid(xl,yl);
surf(X,Y,zeros(size(X)))
shading flat
alpha 0.1
grid on;
grid minor;
title(header2);
xlabel('x in AE')
ylabel('y in AE');
zlabel('z in 0.1 AE');
legend(Planets.Name,'Sonne','Richtung Frühlingspunkt','Venus Transit',...
       'location','northeast');
legend boxoff;
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------