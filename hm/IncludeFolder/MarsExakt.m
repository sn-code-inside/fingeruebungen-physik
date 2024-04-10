% -------------------------------------------------------------------------
% MarsExakt.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial basiert in Teilen 
% auf C++ Programmstrukturen und der Beschreibung/Erläuterung in  
%
% "Astronomie mit dem Personal Computer"
%
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Berechnet die heliozentrischen Koordinaten des Venus.
% Alle Werte beziehen sich auf  Äquinoktium des Datums
%
% MarsPos = MarsExakt(t,MarsData)
%
%  Eingabe:
%    t          Zeit in Julianischen Jahrhunderten seit JD2000
%    MarsData    Eingelesene Tabelle Venus
%  Ausgabe:
%   VenPos.ekl  Heliozentrisch-ekliptikale Koordinaten Venus
%   VenPos.xyz  xyz Koordinaten (heliozentr.) Venus 
% 
% 
% 
% -------------------------------------------------------------------------

function MarsPos = MarsExakt(t,MarsData)

  pi2     = 2*pi;
  Arcs    = 3600*180/pi;
  dl = 0; dr = 0; db = 0;
 
  % Mittlere Anomalien der Planeten in [rad]
  M2 = pi2 * Frac ( 0.1382208 + 162.5482542*t );
  M3 = pi2 * Frac ( 0.9926208 +  99.9970236*t );
  M4 = pi2 * Frac ( 0.0538553 +  53.1662736*t );
  M5 = pi2 * Frac ( 0.0548944 +   8.4290611*t );
  M6 = pi2 * Frac ( 0.8811167 +   3.3935250*t );


  % Stoerungen durch Venus (M4, M2)
  for k =1:12
      i4 = MarsData(k,2);
      i2 = MarsData(k,3);
      ct  = cos(i4*M4+i2*M2); st=sin(i4*M4+i2*M2);
      dlt = MarsData(k,5)*ct+MarsData(k,6)*st;
      drt = MarsData(k,7)*ct+MarsData(k,8)*st;
      dbt = MarsData(k,9)*ct+MarsData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
  end
  % Störungen durch Erde (M4, M3)
  for k =14:66
      i3 = MarsData(k,3);
      i4 = MarsData(k,2);
      iT = MarsData(k,4);
      ct  = cos(i4*M4+i3*M3); st=sin(i4*M4+i3*M3);
      if ~(iT == 0) 
        ct =ct.*t.^iT; st = st.*t.^iT;  
      end
      dlt = MarsData(k,5)*ct+MarsData(k,6)*st;
      drt = MarsData(k,7)*ct+MarsData(k,8)*st;
      dbt = MarsData(k,9)*ct+MarsData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
  end
  % Störungen durch Jupiter (M4, M5)
  for k =68:96
      i4 = MarsData(k,2);
      i5 = MarsData(k,3);
      ct  = cos(i4*M4+i5*M5); st=sin(i4*M4+i5*M5);
      dlt = MarsData(k,5)*ct+MarsData(k,6)*st;
      drt = MarsData(k,7)*ct+MarsData(k,8)*st;
      dbt = MarsData(k,9)*ct+MarsData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
  end
  % Störungen durch Saturn (M4, M6)
  for k =98:109
      i4 = MarsData(k,2);
      i6 = MarsData(k,3);
      ct  = cos(i4*M4+i6*M6); st=sin(i4*M4+i6*M6);
      dlt = MarsData(k,5)*ct+MarsData(k,6)*st;
      drt = MarsData(k,7)*ct+MarsData(k,8)*st;
      dbt = MarsData(k,9)*ct+MarsData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
  end
  % Ekliptikale Koordinaten ([rad],[AE])
  dl  =  dl + 52.49*sin(pi2*(0.1868+0.0549*t)) + ...
         0.61*sin(pi2*(0.9220+0.3307*t)) + ...
         0.32*sin(pi2*(0.4731+2.1485*t)) + ...
         0.28*sin(pi2*(0.9467+0.1133*t));
  dl  =  dl + 0.14 + 0.87*t - 0.11*t.*t;

  % Ekliptikale Koordinaten ([rad],[AE])
  l = pi2*Frac(0.9334591 + M4/pi2 + ( (6615.5+1.1.*t).*t + dl )/1296.0e3);
  r = 1.5303352 + 0.0000131*t  +  dr * 1.0e-6;
  b = ( 596.32 + (-2.92 - 0.10.*t).*t  +  db )/Arcs;

  MarsPos.ekl = [r; l; b];
  MarsPos.xyz = CalcXYZfromAngles(MarsPos.ekl); % xyz helio-ekliptikal
  MarsPos.Name ="Mars";
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------