% -------------------------------------------------------------------------
% PlanetOrbits.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% PlanetOrbits
% Berechnet die x-, y- und z-Koordinaten der Planetenbahnen
% -------------------------------------------------------------------------

function Planeten=PlanetOrbits(BaPa)
% Eingabe: 
%   BaPa    = Bahnparameter (siehe OrbitParameterEinlesen.m und
%               Bahnparameter03.csv, Bahnparameter04.csv und Tabelle 3.8 
% Ausgabe: 
% Planeten  = Planeten mit Name, x-, y- und z-Koordinaten

    % Winkelarray u 
    u=linspace(0,360,361);
    % Dimensionen der Bahnparameter-Tabelle 
    SizeT = size(BaPa);
    
    % Iteration durch die Planeten 
    for iPlot=1:SizeT(1)
        % Exzentrizitaet 
        ePn=BaPa.eP(iPlot);
        % Laenge der großen Halbachse 
        aPn=BaPa.aP(iPlot);
        
        Planets(iPlot).Name=string(BaPa.Name(iPlot));
        % Berechnung der Ellipse
        fac=aPn*(1-ePn*ePn)./(1+ePn*cosd(u));
        x = cosd(u).*fac;
        y = sind(u).*fac;
        z = cosd(u).*0;
        yin=[x;y;z];
        
        % Inklination/Neigung der Bahnebene im Bogenmass
        iPn=deg2rad(BaPa.iP(iPlot));
        % Winkeldistanz vom Fruehlingspunkt bis zum aufsteigenden Knoten 
        % in der Ekliptik in Bogenmass
        OmegaPn=deg2rad(BaPa.OmegaP(iPlot));
        % Perihelargument in Bogenmass 
        PiqPn=deg2rad(BaPa.PiqP(iPlot));
        % Laenge des Perihels in Bogenmass (vgl. Formel 3.66) 
        omegaPn=PiqPn-OmegaPn;
        % Matrixmultiplikation 
        yout=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
        for m=1:3 Planets(iPlot).xyz(m,:)=yout(m,:); end
    end
    Planeten = Planets;
end 
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
