% -------------------------------------------------------------------------
% ZGL02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die ZGL auf Basis einer Reihenentwicklung.
% (Äquinoktium und Ekliptik des Datum J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
dt1 = datetime('2000-01-01 12:00:00');
T = juliandate(dt1);
MJuDa = juliandate(dt1,'modifiedjuliandate');
DatumStr =  string(dt1,'dd.MM.yyyy');
T0=(T-2451545)/36525;
eps=deg2rad(23.43929111-(46.8150+0.00059*T0-0.001813*T0*T0)*T0/3600);
pi2 = 2*pi;
U = 365.25;

%-------------------------------------------------------------------------
%Begin Rechnung

T_vector=(T):(T+365);
MJuDa_vector = (MJuDa):(MJuDa+365);
Timey = 1:366;

% ZGL nach Formel In PhysicsUncut
tauw1(1,:) = -7.36353*sin(0.017202*(Timey)-0.062660);
tauw1(2,:) = -9.93481*sin(0.034404*(Timey)+0.361835);
tauw1(3,:) = -0.329587*sin(0.051606*(Timey)+3.463816-pi);
tauw1(4,:) = -0.21222*sin(0.068808*(Timey)+3.872289-pi);
tauw1(5,:)=zeros(1,366);
for i=1:4 tauw1(5,:)=tauw1(5,:)+tauw1(i,:); end


%------------------------------------------------------------------------------
% Graphische Ausgabe
Nulllinie=zeros(366);

titlestr   =  strings([6,25]);
titlestr(1,:)= ' Formel (4-79) ';
titlestr(2,:)= ' Lösung Rempel ';
titlestr(3,:)= ' JPL';
TimeStr=string(hours(12)+minutes(0)+seconds(0),'hh:mm');
header1='Zeitgleichung Empirische Daten und Keplerlösung für UT='+TimeStr; 
figure('Name',header1);
plot(Timey, tauw1(5,:), 'Color', Colors(3,:),'LineStyle','-','LineWidth',2);
hold on;
plot(Timey, tauw1(1,:), 'Color', Colors(2,:),'LineStyle','-','LineWidth',2);
plot(Timey, tauw1(2,:), 'Color', Colors(4,:),'LineStyle','-','LineWidth',2);
plot(Timey, 10*tauw1(3,:), 'Color', Colors(8,:),'LineStyle','-.','LineWidth',2);
plot(Timey, 10*tauw1(4,:), 'Color', Colors(9,:),'LineStyle','-.','LineWidth',2);
plot(Timey, Nulllinie, 'Color', 'k','LineWidth',1);
hold on;
ylim([-20 20]);
xlim([0 370]);
grid on;
header2=strcat('Aequinoktium :', DatumStr);
xlabel('Tag')
ylabel('ZGL (Anteile) in Minuten');
legend(' ZGL', ' Term 1 ({\omega})',' Term 2 ({2\omega})',' 10 x Term 3 ({3\omega})',' 10 x Term 4 ({4\omega})','Location','southeast');
legend boxoff;
set(gca,'Fontsize',18);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------