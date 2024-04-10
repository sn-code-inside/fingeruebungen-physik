% -------------------------------------------------------------------------
% SonneExakt.m
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
% Berechnet die geozentrischen Koordinaten der Sonnne.
% Alle Werte beziehen sich auf  Äquinoktium des Datums
%  Eingabe:
%    t          Zeit in Julianischen Jahrhunderten seit JD2000
%    epsE       Ekiptikschiefe der Erde in rad
%    rMass      Einheit der Abstandes, entweder 'km', Erdradien 'RE' oder
%               astronomischen Einheiten 'aE'
%    SunData    Eingelesene Tabelle Parameterdaten
%  Ausgabe:
%    SunPos.ekl   Geozentrisch-ekliptikale Koordinaten Sonne
%    SunPos.xyz   xyz Koordinaten (geozentr.) Sonne 
%    SunPos.geo   Geozentrische äquatoriale xyz Koordinaten 
%    SunPos.equ   Geozentrische äquatoriale  Koordinaten in rad
% -------------------------------------------------------------------------

function SunPos = SonneExakt(t,SunData,epsE, rMass)

    pi2     = 2*pi;
    Arcs    = 3600*180/pi;

    dl = 0; dr = 0; db = 0;

    %Mittlere Anomalien
    M2 = pi2 * Frac ( 0.1387306 + 162.5485917*t );
    M3 = pi2 * Frac ( 0.9931266 +  99.9973604*t );
    M4 = pi2 * Frac ( 0.0543250 +  53.1666028*t ); 
    M5 = pi2 * Frac ( 0.0551750 +   8.4293972*t );
    M6 = pi2 * Frac ( 0.8816500 +   3.3938722*t ); 

    D  = pi2 * Frac ( 0.8274 + 1236.8531*t );
    A  = pi2 * Frac ( 0.3749 + 1325.5524*t );      
    U  = pi2 * Frac ( 0.2591 + 1342.2278*t );


    % Störungen durch Venus
    for k =1:23
      iT = SunData(k,4);
      i5 = SunData(k,2);
      i6 = SunData(k,3);
      ct  = cos(i5*M3+i6*M2); st=sin(i5*M3+i6*M2);
      if ~(iT == 0) 
        ct =ct.*t.^iT; st = st.*t.^iT;  
      end
      dlt = SunData(k,5)*ct+SunData(k,6)*st;
      drt = SunData(k,7)*ct+SunData(k,8)*st;
      dbt = SunData(k,9)*ct+SunData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
    end
    % Störungen durch Mars
    for k =25:37
      i5 = SunData(k,2);
      i7 = SunData(k,3);
      ct  = cos(i5*M3+i6*M4); st=sin(i5*M3+i7*M4);
      dlt = SunData(k,5)*ct+SunData(k,6)*st;
      drt = SunData(k,7)*ct+SunData(k,8)*st;
      dbt = SunData(k,9)*ct+SunData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
    end
    % Störungen durch Jupiter
    for k =39:50
      i5 = SunData(k,2);
      i7 = SunData(k,3);
      ct  = cos(i5*M3+i6*M5); st=sin(i5*M3+i7*M5);
      dlt = SunData(k,5)*ct+SunData(k,6)*st;
      drt = SunData(k,7)*ct+SunData(k,8)*st;
      dbt = SunData(k,9)*ct+SunData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
    end
    % Störungen durch Saturn
    for k =52:55
      i5 = SunData(k,2);
      i7 = SunData(k,3);
      ct  = cos(i5*M3+i6*M6); st=sin(i5*M3+i7*M6);
      dlt = SunData(k,5)*ct+SunData(k,6)*st;
      drt = SunData(k,7)*ct+SunData(k,8)*st;
      dbt = SunData(k,9)*ct+SunData(k,10)*st;
      dl  = dl + dlt;
      dr  = dr + drt;
      db  = db + dbt; 
    end

    % Abstand des Schwerpunkts des Erde-Mond-Systems vom Erdmittelpunkt
    dl = dl +  6.45*sin(D) - 0.42*sin(D-A) + 0.18*sin(D+A) + ...
               0.17*sin(D-M3) - 0.06*sin(D+M3);
    dr = dr + 30.76*cos(D) - 3.06*cos(D-A) + 0.85*cos(D+A) - ...
               0.58*cos(D+M3) + 0.57*cos(D-M3);
    db = db + 0.576*sin(U);

    % Langperiodische Stoerungen
    dl = dl  + 6.40 * sin ( pi2*(0.6983 + 0.0561.*t) ) + ...
             1.87 * sin ( pi2*(0.5764 + 0.4174.*t) ) + ...
             0.27 * sin ( pi2*(0.4189 + 0.3306.*t) ) + ...
             0.20 * sin ( pi2*(0.3581 + 2.4814.*t) );
    % Ekliptikale Koordinaten ([rad],[AE])
    l = pi2 * Frac ( 0.7859453 + M3/pi2 + ...
                 ( (6191.2+1.1.*t).*t + dl ) / 1296.0e3 );
    r = 1.0001398 - 0.0000007.* t + dr * 1.0e-6;
    b = db / Arcs;
    % Hier Fallunterscheidung ob in AE oder in Erdradien
    if strcmpi(rMass,'km') 
       fac =149.5e6;  % r in km
    else
      if strcmpi(rMass,'RE') 
         fac = 23454.78; % r in Erdradien
      else
         fac =1;  % r in AE Einheiten
      end
    end
    SunPos.ekl = [r*fac; l ; b]; 
    SunPos.xyz = CalcXYZfromAngles(SunPos.ekl); % xyz helio-ekliptikal
    SunPos.geo = mtimes(R_x(-epsE),SunPos.xyz); % xyz in Equ-System
    SunPos.equ = CalcAnglesfromXYZ(SunPos.geo); % RA, DEC 
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------