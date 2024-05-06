% -------------------------------------------------------------------------
% R_x.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Drehmatrix R_x 
% Verwendung in PQR.m
%
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_x 
% 
% -------------------------------------------------------------------------

function yout=R_x(phi)
   yout=[1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------