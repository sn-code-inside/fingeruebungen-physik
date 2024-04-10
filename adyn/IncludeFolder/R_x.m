%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% R_x.m
% MATLAB Programm zu Fingerübungen fuer Physik von Michael Kaschke
% Kapitel Himmelsmechanik
% Alle Rechte beim Autor
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%

%-------------------------------------------------------------------------
% R_x
% Drehmatrix R_x 
% Verwendung in PQR.m
%-------------------------------------------------------------------------
%%

function yout=R_x(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_x 

   yout=[1 0 0;0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
end   
%%------------------------------------------------------------------------
% Ende Funktion
%-------------------------------------------------------------------------

