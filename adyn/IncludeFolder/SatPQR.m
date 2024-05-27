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
  
  aP = SatSet.BaPa(1);
  eP = SatSet.BaPa(2);
  iP = SatSet.BaPa(3);
  MP = SatSet.BaPa(4);  %Mittlere Anomalie zum Startzeitpunkt
  OmegaP = SatSet.BaPa(5);
  omegaP = SatSet.BaPa(6);
  nP  = sqrt(GME/aP^3)*86400;  %Umlaufgeschwindigkeit in 1/s
% Mittlere Anomalie
  % Berechnung der exzentrischen Anomalie mit Startwert M(1) = MPn
  for k=1:length(T)
      M(k) = deg2rad(MP) + nP*(T(k)-T(1));
      Ecc(k)=EAnom(M(k),eP);
%       fprintf('\n Ecc(k)     = %12.9f     (Exzentr. Anomalie)', Ecc(k));
  end
% Berechnung der Bahndaten
  cosE=cos(Ecc);
  sinE=sin(Ecc);
  fac=sqrt((1-eP)*(1+eP));
  % Bahnkoordinaten 
  dis     = aP.*(1-eP*cosE);
  Vel     = sqrt(GME*aP)./dis;
  uvw     = [aP.*(cosE-eP);aP*fac.*sinE;0.*Ecc]; 
  vel     = Vel.*[-sinE;fac.*cosE;0.*Ecc];
  iP=deg2rad(iP);
  OmegaP=deg2rad(OmegaP);
  omegaP=deg2rad(omegaP);
  M_PQR = PQR(-OmegaP, -iP,-omegaP); 
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