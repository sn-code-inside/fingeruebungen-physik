% ------------------------------------------------------------------------
% VektorAnalysis01m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Bewegung einer Achterbahn in einem
% kontinuierlichen Modell eines elastischen Bandes.
% -------------------------------------------------------------------------

clc
clear
x = -2*pi:0.5:2*pi; %sets up an x-axis
y = -2*pi:0.5:2*pi; %sets up a y-axis
[x,y] = meshgrid(x,y); %takes the axis info and makes a grid
%of x and y values
u = -sin(y); %calculates the u velocity
v = sin(x); %calculates the v velocity
quiver(x,y,u,v) %draws the velocity vectors on the grid
xlabel('x-axis')
ylabel('y-axis')
title('velocity vector plot')