% -------------------------------------------------------------------------
% SatPQR.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Astrodynamik aus
% "Finger√ºbungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die geozentrisch-‰quatorialen Koordinaten von Satelliten aus
% den Bahnpaarmetern f√ºr eine bestimmte Zeit 
%
% Eingabe:
%   T                               = Zeit in Julianischem Datum
%   aPn,ePn,iPn,MPn,OmegaPn,omegaPn = Orbitparameter in km
%                                     bzw. ∞ √ºbergeben
%  
% Ausgabe:
%  xyz      = geozentrisch-‰quatorialen Koordinaten (xyz) des Satelliten
%  vel      = Gesschwindigkeitskomponenten (xyz) des Satelliten
%  Sat.az; Sat.el, Sat.r   = geozentrisch-‰quatorialen Koordinaten
%                            (alpha, delta, r) des Satelliten
%  Time     = Zeitst√ºtzstellen in JD
%  Name     = Name des Satelliten
% -------------------------------------------------------------------------

function SatData= SatPQR(T, SatSet, GME) 
  % Zeit  T in Julianischen Datum 
  
  aPn = SatSet.BaPa(1);
  ePn = SatSet.BaPa(2);
  iPn = SatSet.BaPa(3);
  MPn = SatSet.BaPa(4);  %Mittlere Anomalie zum Startzeitpunkt
  OmegaPn = SatSet.BaPa(5);
  omegaPn = SatSet.BaPa(6);
 
  nP  = sqrt(GME/aPn^3)*86400;  %Umlaufgeschwindigkeit in 1/s
% Mittlere Anomalie
  M = deg2rad(MPn) + nP*T;
  % Berechnung der exzentrischen Anomalie
  Ecc=EAnom(M,ePn);
% Berechnung der Bahndaten
  cosE=cos(Ecc);
  sinE=sin(Ecc);
  fac=sqrt((1-ePn)*(1+ePn));
  % Bahnkoordinaten 
  dis     = aPn.*(1-ePn*cosE);
  Vel     = sqrt(GME*aPn)./dis;
  uvw     = [aPn.*(cosE-ePn);aPn*fac.*sinE;0.*Ecc]; 
  vel     = Vel.*[sinE;fac.*cosE;0.*Ecc];
  iPn=deg2rad(iPn);
  OmegaPn=deg2rad(OmegaPn);
  omegaPn=deg2rad(omegaPn);
  M_PQR = PQR(-OmegaPn, -iPn,-omegaPn); 
  Sat.xyz = mtimes(M_PQR,uvw);
  Sat.vel = mtimes(M_PQR,vel);
  xP = Sat.xyz(1,:);
  yP = Sat.xyz(2,:);
  zP = Sat.xyz(3,:);
 [Sat.az,Sat.el,Sat.r] = cart2sph(xP,yP,zP);
  Sat.Time = T; 
  Sat.Name = SatSet.Name;
  SatData=Sat;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------