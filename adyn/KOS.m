% -------------------------------------------------------------------------
% KOS.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Koordinaten-Transformation
% 
% -------------------------------------------------------------------------
clc
clear all
close all 
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Definiert symbolischen Koordinaten 
syms r phi lam rvec S  'real'

% Drehung um x-Achse 
rvec   = r*[cos(phi)*cos(lam);cos(phi)*sin(lam);sin(phi)]

e_r    = diff(rvec,r)
e_phi  = (1/r)*diff(rvec,phi)
e_lam  = (1/r/cos(phi))*diff(rvec,lam)


% Koordinatentarnsformations-Matrix

S = [e_r, e_phi, e_lam]


