% -------------------------------------------------------------------------
% R_y.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Drehmatrix R_y 
% Verwendung in PQR.m
% -------------------------------------------------------------------------

function yout=R_y(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_y 

   yout=[cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
