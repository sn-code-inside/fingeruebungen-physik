% -------------------------------------------------------------------------
% R_x.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Drehmatrix R_x 
% Verwendung in PQR.m
% -------------------------------------------------------------------------

function yout=R_x(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_x 

   yout=[1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
end   
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------

