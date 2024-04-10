%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% R_z.m
% MATLAB Programm zu Fingerübungen fuer Physik von Michael Kaschke
% Kapitel Himmelsmechanik
% Alle Rechte beim Autor
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%

%-------------------------------------------------------------------------
% R_z
% Drehmatrix R_z 
% Verwendung in PQR.m
%-------------------------------------------------------------------------
%%

function yout=R_z(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_z 
   yout=[cos(phi) sin(phi) 0 ;...
        -sin(phi) cos(phi) 0 ;...
        0 0 1];
end
%%------------------------------------------------------------------------
% Ende Funktion
%-------------------------------------------------------------------------