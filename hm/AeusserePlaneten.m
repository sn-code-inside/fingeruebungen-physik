% -------------------------------------------------------------------------
% AuesserePlaneten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die Bahnen der äußeren Planeten auf Basis der Keplergleichung 
% für den Zeitraum 1920 -- 2020. 
% (Aequinoktium und Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
dt1 = datetime('1920-01-01 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

dt2 = datetime('2020-01-01 00:00:00');
T2 = juliandate(dt2); % Julianisches Datum  
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

% Einlesen der Bahnparameter, 
% BaPadot bezeichnet die erste zeitliche Ableitung 
dt0 = datetime('2000-01-01 12:00:00');
T0 = juliandate(dt0); % Julianisches Datum
Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T0,Aequi);

%-------------------------------------------------------------------------
% Beginn Rechnung
NP =2000; %Stuetzstellen
T_vector=linspace(T1,T2,NP);
% Berechnung nach Keplerlösung
for iPlot=5:9
    Planets(iPlot)=PlanetPQR(T_vector, BaPa, BaPadot, iPlot);
end

%-------------------------------------------------------------------------
% Graphische Ausgabe

header1='Bahnen der äusseren Planeten';
figure('Name',header1);
subplot(1,2,1);
% Planetenbahnen
for iPlot = 5:9
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', ...
         Colors(iPlot,:),'LineWidth',1);
    hold on
end
axis equal
SonneFP;
for iPlot = 5:9
    % Endpunkt T2
    plot(Planets(iPlot).xyz(1,NP),Planets(iPlot).xyz(2,NP),'d','Color', ...
         Colors(iPlot,:),'LineWidth',2);
    % Startpunkt T1
    plot(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(2,1),'+','Color', ...
         Colors(iPlot,:),'LineWidth',2);
    axis equal
end
ylim([-40 40])
xlim([-40 40]);
grid on;
grid minor;
header2 = strjoin(['Projektion auf die Ekliptik: ',DatumStr1,' - ',...
                  DatumStr2]);
title(header2,'FontSize',14);
xlabel('x in AE')
ylabel('y in AE')
legend(Planets(5).Name,Planets(6).Name,Planets(7).Name,Planets(8).Name,...
       Planets(9).Name,'Sonne',...
       'Richtung Äquinoktium','location','northeast');
legend boxoff;
set(gca,'FontSize',14);

%------------------------------------------------------------------------------

subplot(1,2,2);
% Planetenbahnen
for iPlot = 5:9
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(3,:)*3,'Color', ...#
         Colors(iPlot,:),'LineWidth',1);
    hold on
end
axis equal
SonneFP;
for iPlot = 5:9
   % Endpunkt T2 
   plot(Planets(iPlot).xyz(1,NP),Planets(iPlot).xyz(3,NP)*3,'d','Color',...
        Colors(iPlot,:),'LineWidth',2);
   % Startpunkt T1
   plot(Planets(iPlot).xyz(1,1),Planets(iPlot).xyz(3,1)*3,'+','Color',...
        Colors(iPlot,:),'LineWidth',2);
end
xlim([-40 40]);
ylim([-40 40]);
grid on;
grid minor;
header2 = strjoin(['Senkrecht zur Ekliptik: ',DatumStr1,' - ',...
                  DatumStr2,'']);
title(header2);
xlabel('x in AE')
ylabel('3 z in AE (gestreckt)');
legend(Planets(5).Name,Planets(6).Name,Planets(7).Name,Planets(8).Name,...
       Planets(9).Name,'Sonne','Ebene der Ekliptik','location','northeast');
legend boxoff;
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------