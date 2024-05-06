% -------------------------------------------------------------------------
% Convert2Equ.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Umrechnung von heliozentrischen in geozentrische Koordinaten und weitere
% Umrechnung in aequatoriale Koordinaten 
% -------------------------------------------------------------------------

function PlanetDataOut = Convert2Equ(PlanetIn, Erde, epsrad)
% Eingabe: 
%   PlanetIn    	= Bahnparameter des Planeten
%   Erde            = Bahnparameter der Erde
%   epsrad          = Ekliptikschiefe [eps] = rad 
% Ausgabe: 
%   PlanetDataOut   = Bahnparameter des Planeten inkl. den aequatorialen
%                     Koordinaten 
    if string(PlanetIn.Name) == 'Erde'
        PlanetIn.geo  = -Erde.xyz;    % geozentrisch Sonne, kartesisch
    else
        PlanetIn.geo  = -Erde.xyz + PlanetIn.xyz; % geozentrisch,kartesisch
    end
    PlanetIn.equ=rad2deg(CalcAnglesfromXYZ((mtimes(R_x(-eps),...
                             PlanetIn.geo)))); % geozentrisch,aequatorial 
    PlanetDataOut = PlanetIn;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
