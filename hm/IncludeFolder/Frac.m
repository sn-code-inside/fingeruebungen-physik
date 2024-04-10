% -------------------------------------------------------------------------
% Frac.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Nachkommastellen 
% -------------------------------------------------------------------------

function y = Frac(x)          
% Eingabe: 
%   x   = Dezimalzahl z.B. Ergebnis einer Division
% Ausgabe: 
%   y   = Nachkommastellen 
      y = x-fix(x);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------