% -------------------------------------------------------------------------
% DD.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Ermittelt Grad in Dezimalform (°) aus Grad, Minuten und Sekunden 
% -------------------------------------------------------------------------

function xDD = DD(DMS)
% Eingabe:
%   DG        Grad in dezimaler Darstellung
% Ausgabe:  x[G,M,S]
%   G         Grad
%   M         Minuten
%   S         Sekunden 
  VZ = sign(DMS(1));
  xDD = VZ*DMS(1)+DMS(2)/60 + DMS(3)/3600;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
