% -------------------------------------------------------------------------
% PlanetExzentrizitaeten.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Darstellung der Exzentrizitäten der Planetenbahnen (auf große 
% Halbachse a normiert).
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Einlesen der Bahnparameter, 
% BaPadot bezeichnet die erste zeitliche Ableitung 
T0 = 2451545.0;
Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T0,Aequi);

exz = BaPa.eP; 
u = linspace(0,360,361);

% cosd, sind = Cosinus, Sinus in Grad (degree)
cosu=cosd(u);
sinu=sind(u);
% Iteration durch alle Planeten 
gamma = 1-exz.*exz;

for m=1:9    
    x(m,:)=gamma(m)*(cosu./(1+exz(m).*cosu));
    y(m,:)=gamma(m)*(sinu./(1+exz(m).*cosu));  
end

%------------------------------------------------------------------------------
% Graphische Ausgabe

header1='Keplerlösung'; 
figure('Name',header1);
hold on
for m=1:9
    plot(x(m,:),y(m,:),'Color',Colors(m,:));
end
plot(0,0,'+','Color',Colors(10,:),'LineWidth',5);
ylim([-1.25 1.25]);
xlim([-1.35 1.15]);
axis square;
title('\it r/a \rm Verhältnis der Planetenbahnen');
grid on;
legend(BaPa.Name(:),'location','bestoutside');
legend('boxoff');
axis equal;

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
