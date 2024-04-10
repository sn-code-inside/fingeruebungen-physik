% -------------------------------------------------------------------------
% PlanetPQR.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die heliozentrisch-ekliptikalen Koordinaten der Planeten aus
% den Bahnpaarmetern für eine bestimmte Zeit 
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

function PlanetData= PlanetPQR(T, BaPa, BaPadot, k) 
  pi2=2*pi;
  t=Jd2JJht(T); % Zeit in Julianischen Jahrhunderten seit JD2000
  
% Planetenparameter
  aPn=BaPa.aP(k);
  ePn=BaPa.eP(k);
  NPn=BaPa.NP(k);
  
% Mittlere Anomalie    
  M = deg2rad((BaPa.L0P(k)-BaPa.PiqP(k)) + NPn*t);
  
% Berechnung der exzentrischen Anomalie
  Ecc=EAnom(M,ePn);
% Berechnung der Bahndaten
  cosE=cos(Ecc);
  sinE=sin(Ecc);
  fac=sqrt((1-ePn)*(1+ePn));
  y=aPn*fac*sinE; 
  x=aPn.*(cosE-ePn);
  z=0*x;
  % Bahnkoordinaten 
  uvw  = [x;y;z]; 

  iPn=deg2rad(BaPa.iP(k));
  PiqPn=deg2rad(BaPa.PiqP(k));
  OmegaPn=deg2rad(BaPa.OmegaP(k));
  omegaPn=PiqPn-OmegaPn;
    
  Planet(k).xyz = mtimes(PQR(-OmegaPn,-iPn,-omegaPn),uvw);
  Planet(k).ekl = CalcAnglesfromXYZ(Planet(k).xyz);
  Planet(k).Time = T; 
  Planet(k).Name = BaPa.Name{k};
  PlanetData=Planet(k);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------