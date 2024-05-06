% -------------------------------------------------------------------------
% EpsErde.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Achsneigung der Erde
% -------------------------------------------------------------------------

function epsilon =EpsErde(T0)
% Eingabe: 
%   T0          = Zeit als Julianisches Datum 
% Ausgabe: 
%   epsilon     = Achsneigung in Grad
   T0=Jd2JJht(T0);
   epsilon =23.43929111-(46.8150+0.00059*T0-0.001813*T0*T0)*T0/3600;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------