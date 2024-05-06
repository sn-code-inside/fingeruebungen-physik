% -------------------------------------------------------------------------
% JJht2Jd.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bestimmt das Julianische Datum aus dem Modifiziertem JuDatum
% -------------------------------------------------------------------------

function yout=Jd2JJht(T)
  yout=(T-2451545.0)/36525;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------