% -------------------------------------------------------------------------
% CalcXYZfromAngles.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die kartesischen Koordinaten aus sphaerischen Koordinaten
% -------------------------------------------------------------------------

function VecOut=CalcXYZfromAngles(VecIn) 
% Eingabe: 
%   VecIn   :  sphaerische Koordinaten als Vektor [r theta phi]
% Ausgabe:
%   VecOut  :  kartesische Koordinaten als Vektor [X Y Z]
    VecOut=[VecIn(1,:).*cos(VecIn(2,:)).*cos(VecIn(3,:));
            VecIn(1,:).*sin(VecIn(2,:)).*cos(VecIn(3,:));
            VecIn(1,:).*sin(VecIn(3,:))];
end

%-------------------------------------------------------------------------
% Ende Funktion
%-------------------------------------------------------------------------
