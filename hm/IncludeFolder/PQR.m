% -------------------------------------------------------------------------
% PQR.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Multiplikation mit Transformationsmatrix aus Gauss-Vektoren
% vgl. Kapitel 2.5.3.1 Berechnung der Planetenbahnen mit Gausschen Vektoren
%
% Eingabe: 
%   Omega     = Knotenlage 
%   i         = Bahnneigung zur Ekliptik 
%   w         = Perihellage
% Ausgabe: 
%   yout      = Matrixprodukt 
% -------------------------------------------------------------------------

function yout=PQR(Omega,i,w)

   ytemp=mtimes(R_x(i),R_z(w));
   yout=mtimes(R_z(Omega),ytemp);
   
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------