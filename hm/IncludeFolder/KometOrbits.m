% -------------------------------------------------------------------------
% KometOrbits.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Asteroiden- und Kometenbahn (Ellipsengleichung)
% 
% Eingabe:
%     BaPaK     = Bahnparameter der Kometen 
% Ausgabe: 
%     Kometen   = Kometen mit Name, x-, y-, und z-Koordinate
% -------------------------------------------------------------------------

function Kometen=KometOrbits(BaPaK)
    % Winkelarray u 
    u=linspace(0,360,361);
    % Dimension der Bahnparameter-Tabelle
    SizeT = size(BaPaK);
    
    % Iteration durch die Kometen 
    for iPlot=1:SizeT(1)
        % Exzentrizitaet 
        ePn=BaPaK.eP(iPlot);
        Komets(iPlot).Name=string(BaPaK.Name(iPlot)); 
        % Periheldistanz 
        qPn=BaPaK.qP(iPlot);
        % Berechnung der Ellipse 
        fac=qPn*(1+ePn)./(1+ePn*cosd(u));
        x = cosd(u).*fac;
        y = sind(u).*fac;
        z = cosd(u).*0;
        yin=[x;y;z];
        % Inklination/Neigung der Bahnebene im Bogenmass
        iPn=deg2rad(BaPaK.iP(iPlot));
        % Winkeldistanz vom Fruehlingspunkt bis zum aufsteigenden Knoten 
        % in der Ekliptik in Bogenmass 
        OmegaPn=deg2rad(BaPaK.OmegaP(iPlot));
        % Perihelargument in Bogenmass
        omegaPn=deg2rad(BaPaK.omegaP(iPlot));
        % PQR-Matrixmultiplikation 
        yout=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
        for m=1:3 
            Komets(iPlot).xyz(m,:) = yout(m,:); 
        end
    end
    Kometen=Komets;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------