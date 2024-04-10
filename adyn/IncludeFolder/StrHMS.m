%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%
% StrHMS.m
% MATLAB Programm zu Fingeruebungen fuer Physik von Michael Kaschke
% Kapitel Himmelsmechanik
% Alle Rechte beim Autor
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% Das vorliegende Programmaterial ist mit keiner Verpflichtung oder 
% Garantie irgendwelcher Art verbunden. Die Benutzung der Software ist nur 
% zu privaten oder nicht-kommerziellen Zwecken, wie Bildung und Lehre 
% gestattet.
%%-----------------------------------------------------------------------%%
%%-----------------------------------------------------------------------%%

%-------------------------------------------------------------------------
% StrHMS
% Umrechnung eines Winkels von dezimaler Darstellung zu Stunde, Bogenminute 
% und Bogensekunde
%--------------------------------------------------------------------------
%%

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
%%-------------------------------------------------------------------------
% Ende Funktion
%--------------------------------------------------------------------------

