% -------------------------------------------------------------------------
% Stumpff.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Stumpffunktionen zur Berechnung von Parabel-Naeherungsloesungen der 
% Kometenbahnen
%
% Eingabe: 
%   EA          = Quadrat der exzentrischen Anomalie 
% Ausgabe: 
% c1, c2, c3    = Ergebnis der Stumpfffunktionen 1 bis 3 
% -------------------------------------------------------------------------

function [c1,c2,c3]=Stumpff(EA) 

% Initialisierung 
  cx1=0;
  cx2=0;
  cx3=0;
  add=1;
  n=1;
% Approximation durch Iteration (vgl. Formel 3.137) 
  for n=1:100
      cx1=cx1+add;
      add=add/(2*n);
      cx2=cx2+add;
      add=add/(2*n+1);
      cx3=cx3+add;
      add=add.*(-EA);
  end
% Ausgabe 
  c1=cx1;
  c2=cx2;
  c3=cx3;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------