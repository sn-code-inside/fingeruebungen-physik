% -------------------------------------------------------------------------
% SonnePQR.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die heliozentrisch-ekliptikalen Koordinaten einens Quasiplaneten 
% aus den Bahnparametern, wenn dieser sich auf einer
% Kreisbahn entlang seiner Aequtorebene um die Sonnen bewegen wuerde.
% Diese Koordinaten braucht man in der Uebung, wenn man die RA der
% mittleren Sonne von einem Planeten aus berechnen will.
%
% Eingabe:
%   T                 = Zeit in Julianischem Datum
%   BaPa, BaPadot     = Orbitparameter im Grad übergeben
%   k                 = Planetennummer
% Ausgabe:
%  ekl      = Heliozentrisch-ekliptikale Koordinaten des Planeten k
%  xyz      = Heliozentrisch-ekliptikale Koordinaten (xyz) des Planeten k
%  Time     = Zeitstützstellen in JD
%  Name     = Name des Planeten
% -------------------------------------------------------------------------

function PlanetData= SonnePQR(T, BaPa, BaPadot, k) 

  pi2=2*pi;
  t=Jd2JJht(T); % Zeit in Julianischen Jahrhunderten seit JD2000
% Planetenparameter
  aPn=BaPa.aP(k);
  L0Pn=BaPa.L0P(k);
  M0Pn=BaPa.M0P(k);
  NPn=BaPa.NP(k);
  ePn=0;  % Einfach Exzentrizität gleich 0 gesetzt
% Mittlere Anomalie - Variante 1 der Berechnung
  M  = pi2 * Frac((M0Pn + NPn*t)/360);
  y=aPn*sin(M); 
  x=aPn*cos(M);
  z=0*x;
% Bahnkoordinaten 
  uvw  = [x;y;z]; 
  %   Hier Addition der Bahnneigungen und Achsneigung, bei der Erde in 
  %   entgegengesetzter Richtung, 
  %   Eventuell müssen wir auch noch eine Drehung berücksichtigen.
  if k == 3 
      iPn=deg2rad(BaPa.iP(k)) - deg2rad(EpsErde(T));
  else
      iPn=deg2rad(BaPa.iP(k)) + deg2rad(BaPa.epsP(k));
  end
  PiqPn=deg2rad(BaPa.PiqP(k));
  OmegaPn=deg2rad(BaPa.OmegaP(k));
  omegaPn=PiqPn-OmegaPn;
%   Planet(k).ekl  = [aPn; M; z];
%   Planet(k).xyz  = CalcXYZfromAngles(Planet(k).ekl);
  Planet(k).xyz  = mtimes(PQR(-OmegaPn,-iPn,-omegaPn),uvw);
  Planet(k).ekl  = CalcAnglesfromXYZ(Planet(k).xyz);
  Planet(k).Time = T; 
  Planet(k).Name = BaPa.Name{k};
  PlanetData=Planet(k);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------