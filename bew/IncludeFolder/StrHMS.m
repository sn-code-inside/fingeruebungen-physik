% -------------------------------------------------------------------------
% StrHMS.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Umrechnung eines Winkels von dezimaler Darstellung zu Stunde, Bogenminute 
% und Bogensekunde
% -------------------------------------------------------------------------

function xStrHMS = StrHMS(DH)
% Eingabe:
%   DH        Stunde in dezimaler Darstellung
%
% Ausgabe:  xHMS[H,M,S]
%   H         Stunde
%   M         Minuten
%   S         Sekunden
  if isnan(DH)
    xStrHMS = '--------';
  else
    Sec  = 3600*DH;
    xHMS = [floor(Sec/3600),floor(rem(Sec,3600)/60),rem(rem(Sec,3600),60)];
    xStrHMS = string(duration(xHMS),'hh:mm:ss.SSS'); 
  end
end
% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
