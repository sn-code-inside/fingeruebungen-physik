% -------------------------------------------------------------------------
% Refraktionsformeln.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnete Normalrefraktion.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

h=linspace(-10,90,1000);
R1=1.02./tand(h+10.3./(h+5.11));
R2=1.0./tand(h+7.31./(h+4.4));
Nulllinie=zeros(1000);

figure('Name','Refraktion');
plot(h,R1,':',h,R2,'LineWidth',2);
hold on;
plot(h,Nulllinie,'k','LineStyle',':','LineWidth',2);
title('Normalrefraktion');
grid on;
grid minor;
xlim([0 10]);
ylim([0 40]);
grid on;
xlabel('\it h, h´ \rm in °')
ylabel('\it \color {blue} R(h), \color {red} R(h´) \rm \color {black} in Bogenminuten');
legend('\it R(h)', '\it R(h´)');
legend boxoff;
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------