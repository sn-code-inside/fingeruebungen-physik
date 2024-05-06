% -------------------------------------------------------------------------
% R_z.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Drehmatrix R_z 
% Verwendung in PQR.m
%
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_z 
% 
% -------------------------------------------------------------------------

function yout=R_z(phi)
   yout=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
