%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% R_y.m
% MATLAB Programm zu Fingerübungen fuer Physik von Michael Kaschke
% Kapitel Himmelsmechanik
% Alle Rechte beim Autor
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%

%-------------------------------------------------------------------------
% R_y
% Drehmatrix R_y 
%-------------------------------------------------------------------------
%%

function yout=R_y(phi)
% Eingabe:  
%   phi     = Rotationswinkel im Bogenmaß [rad]
% Ausgabe:
%   yout    = Rotationsmatrix R_y 

   yout=[cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
end

%%------------------------------------------------------------------------
% Ende Funktion
%-------------------------------------------------------------------------
