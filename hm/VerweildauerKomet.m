% -------------------------------------------------------------------------
% VerweildauerKomet.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Verweildauer eines Kometen in Erdnähe.
% -------------------------------------------------------------------------

clc
clear all
close all 

k=linspace(0,1,100);

T = 4*365.25*sqrt(1+3*k-4*k.^3)/(3*sqrt(2)*2*pi);


figure('Name','Verweildauer');
plot(k,T,'LineWidth',2);
hold on;
grid on;
xlim([0 1]);
ylim([0 100]);
ylabel('\it T_{ges}    in Tagen')
xlabel(' Perihelabstand in AE');
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
