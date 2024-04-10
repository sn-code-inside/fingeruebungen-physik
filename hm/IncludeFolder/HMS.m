% -------------------------------------------------------------------------
% HMS.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Ermittelt Stunde, Minuten und Sekunden zur dezimalen Stunde 
% -------------------------------------------------------------------------

function xHMS = HMS(DH)
% Eingabe:
%   DH        Stunde in dezimaler Darstellung
%
% Ausgabe:  xHMS[H,M,S]
%   H         Stunde
%   M         Minuten
%   S         Sekunden
  Sec  = 3600*DH;
  xHMS = [floor(Sec/3600),floor(rem(Sec,3600)/60),floor(rem(rem(Sec,3600),60))];
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------