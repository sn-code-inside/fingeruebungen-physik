% -------------------------------------------------------------------------
% PlanetPQR2.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die heliozentrisch-ekliptikalen Koordinaten der Planeten aus
% den Bahnpaarmetern für eine bestimmte Zeit und auch noch die
% Geschwindigkeiten
%
% Eingabe:
%   T                 = Zeit in Julianischem Datum
%   BaPa, BaPadot     = Orbitparameter im Grad übergeben
%   k                 = Planetennummer
%   muG               = G*MS
% Ausgabe:
%  ekl      = Heliozentrisch-ekliptikale Koordinaten des Planeten k
%  vel      = Heliozentrisch-kartesische Geschwindigkeitn des Planeten k
%  xyz      = Heliozentrisch-kartesische Koordinaten (xyz) des Planeten k
%  Time     = Zeitstützstellen in JD
%  Name     = Name des Planeten
% -------------------------------------------------------------------------

function PlanetData= PlanetPQR2(T, BaPa, BaPadot, k, muG) 
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
  dis  = aPn.*(1-ePn*cosE);
  Vel  = sqrt(muG*aPn)./dis;
  vel  = Vel.*[-sinE;fac.*cosE;0.*Ecc];
  
  iPn=deg2rad(BaPa.iP(k));
  PiqPn=deg2rad(BaPa.PiqP(k));
  OmegaPn=deg2rad(BaPa.OmegaP(k));
  omegaPn=PiqPn-OmegaPn;
  M_PQR         = PQR(-OmegaPn,-iPn,-omegaPn);
  Planet(k).xyz = mtimes(M_PQR,uvw);
  Planet(k).vel = mtimes(M_PQR,vel);
  Planet(k).ekl = CalcAnglesfromXYZ(Planet(k).xyz);
  Planet(k).Time = T; 
  Planet(k).Name = BaPa.Name{k};
  PlanetData=Planet(k);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------