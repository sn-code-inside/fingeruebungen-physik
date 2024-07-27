% -------------------------------------------------------------------------
% FassFormelGauss.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Volumen von Rotationskörpern
% über Gaussschen Integralsatz
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')

h = 2;
a = 1.5;
b = h/sqrt(1-1/a^2);

z =linspace(-h,h,201);

rho = a*sqrt(1-z.^2/b^2);

V  = 2*pi*a^2*h*(1-h^2/b^2/3);

figure()
% plot(rho,z)

phi=0:0.1:2*pi+0.1;

[Z,PHI]=meshgrid(z,phi);

surf(rho.*cos(PHI),rho.*sin(PHI),Z);

shading interp;
camlight;
axis equal


