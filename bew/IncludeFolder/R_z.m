% -------------------------------------------------------------------------
% R_z.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Drehmatrix R_z 
% Verwendung in PQR.m
% -------------------------------------------------------------------------

function yout=R_z(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_z 

   yout=[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0; 0 0 1];
end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------