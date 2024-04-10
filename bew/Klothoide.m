% -------------------------------------------------------------------------
% Klothoide.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet eine Klothoide als Kurvenfunktion
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
p      = 1;                         % Klothoidenparameter
NPoints = 1000;
t     = linspace(0,2*pi,NPoints);

% Funktionen
fun1 = @(x) cos((pi/2)*x.^2);
fun2 = @(x) sin((pi/2)*x.^2);

for k=1:NPoints
    x(k) = p*sqrt(pi).*integral(fun1,0,t(k));
    y(k) = p*sqrt(pi).*integral(fun2,0,t(k));
end

[My,kmax] =max(y);
Mx        =x(kmax);

for k=1:kmax
    x2(k) = 2*x(kmax)+p*sqrt(pi).*integral(fun1,0,-t(k));
    y2(k) = -p*sqrt(pi).*integral(fun2,0,-t(k));
end

% Darstellung
figure(1);
plot (x,y);
xlabel('x')
ylabel('y')
title('Klothoide')
grid on


figure(2)
plot (x(1:kmax),y(1:kmax))
hold on;
plot (x2,y2)
xlabel('x')
ylabel('y')
title('Looping aus Klothoidenaesten')
grid on


