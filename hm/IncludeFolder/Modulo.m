% -------------------------------------------------------------------------
% Modulo.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Modulo-Operation 
% -------------------------------------------------------------------------

function z = Modulo(x,y)
% Eingabe: 
%   x   = Dividend
%   Y   = Divisor
% Ausgabe: 
%   z   = Modulo 
 z = y*Frac(x/y);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
