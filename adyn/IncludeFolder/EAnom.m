% -------------------------------------------------------------------------
% EAnom.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die exzentrische Anomalie aus M und e
% -------------------------------------------------------------------------

function Ecc = EAnom(M,e)
    f= 100000;
    if e < 0.8
        E=M;
    else
        E=pi;
    end
    while abs(f) > 1e-11
       f=E-e.*sin(E)-M;
       E=E-f./(1-e.*cos(E));
    end
    Ecc=E;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
