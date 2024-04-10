% -------------------------------------------------------------------------
% MondMeeus.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die geozentrischen Koordinaten des Mondes.
% Basis Brownsche Mondtheorie nach Jean Meeus.
% Alle Werte beziehen sich auf  Äquinoktium des Datums.
% 
%  t     : Zeit in Julianischen Jahrhunderten seit JD2000
%  epsE  : Ekliptikschiefe in Bogenmaß
%  rMass : ANgabe für Masseinheit für r
%  Ausgabe:
%  MoonPos.ekl      Geozentrische Ekliptikale Korordinaten des Mondes
%  MoonPos.xyz      xyz Koordinaten des Mondes bezogen auf Erdmittelpunkt
%  MoonPos.geo      Geozentrische äquatoriale xyz Koordinaten des Mondes
%  MoonPos.equ      Geozentrische äquatoriale Koordinaten des Mondes 
% -------------------------------------------------------------------------

function MoonPos = MondMeeus(t,epsE,rMass)  
  pi2 = 2*pi;
  t2  = t.*t;
  t3  = t.*t2;
  t4  = t.*t3;
  % mittl. Länge Mond
  L0  =  wrapTo360(218.3164477 + 481267.88123421.*t - ...
         0.0015786.*t2 + t3/538841- t4/65194000);  
  % mittl. Abstand Mond_Sonne 
  D   =  wrapTo360(297.8501921 + 445267.1114034.*t - ...
         0.0018819.*t2 + t3/545868 - t4/113065000);  
  % mittl. Anomalie Sonne
  M   =  wrapTo360(357.5291092 +  35999.0502909.*t - ...
         0.0001536.*t2 + t3/24490000);  
  % mittl. Anomalie Mond
  MS  =  wrapTo360(134.9633964 + 477198.8675055.*t + ...
         0.0087414.*t2 - t3/69699 + t4/14712000);  
  % mittl. Abstand vom Knoten   
  F   =  wrapTo360(93.2720950 + 483202.0175233.*t - ...
         0.0036539.*t2 - t3/3526000 + t4/863310000);  
  % Planetare Störungen   
  A1 = wrapTo360(119.75 + 131.849.*t);
  A2 = wrapTo360( 53.09 + 479264.290.*t);
  A3 = wrapTo360(313.45 + 481266.484.*t);
  EP  = 1 - 0.002516.*t - 0.0000074.*t2;
  % Anfnagswerte Störungen
  DLR = [0 0]; DB = 0;
  % Berechne Störungen 
    DLR = ADDdldr(0,0,1,0,6288774,-20905355,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,-1,0,1274027,-3699111,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,0,0,658314,-2955968,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,2,0,213618,-569925,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,1,0,0,-185116,48888,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,0,2,-114332,-3149,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,-2,0,58793,246158,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-1,-1,0,57066,-152138,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,1,0,53322,-170733,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-1,0,0,45758,-204586,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,1,-1,0,-40923,-129620,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,0,0,0,-34720,108743,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,1,1,0,-30383,104755,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,0,-2,15327,10321,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,1,2,-12528,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,1,-2,10980,79661,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,0,-1,0,10675,-34782,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,3,0,10034,-23210,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,0,-2,0,8548,-21636,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,1,-1,0,-7888,24208,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,1,0,0,-6766,30824,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,0,-1,0,-5163,-8379,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,1,0,0,4987,-16675,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-1,1,0,4036,-12831,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,2,0,3994,-10445,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,0,0,0,3861,-11650,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,-3,0,3665,14403,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,1,-2,0,-2689,-7003,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,-1,2,-2602,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-1,-2,0,2390,10056,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,0,1,0,-2348,6322,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-2,0,0,2236,-9884,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,1,2,0,-2120,5751,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,2,0,0,-2069,0, EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-2,-1,0,2048,-4950,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,1,-2,-1773,4130,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,0,2,-1595,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,-1,-1,0,1215,-3958,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,2,2,-1110,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(3,0,-1,0,-892,3258,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,1,1,0,-810,2616,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,-1,-2,0,759,-1897,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,2,-1,0,-713,-2117,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,2,-1,0,-700,2354,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,1,-2,0,691,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-1,0,-2,596,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,0,1,0,549,-1423,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,4,0,537,-1117,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,-1,0,0,520,-1571,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,0,-2,0,-487,-1739,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,1,0,-2,-399,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,0,2,-2,-381,-4421,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,1,1,0,351,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(3,0,-2,0,-340,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(4,0,-3,0,330,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,-1,2,0,327,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(0,2,1,0,-323,1165,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(1,1,-1,0,299,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,3,0,294,0,EP,D,M,MS,F,DLR);
    DLR = ADDdldr(2,0,-1,-2,0,8752,EP,D,M,MS,F,DLR);
    
    DB = ADDdb(0,0,0,1,5128122	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,1,1,280602	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,1,-1,277693	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,0,-1,173237	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,-1,1,55413	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,-1,-1,46271	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,0,1,32573	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,2,1,17198	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,1,-1,9266	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,2,-1,8822	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,0,-1,8216	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,-2,-1,4324	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,1,1,4200	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,1,0,-1,-3359	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,-1,1,2463	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,0,1,2211	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,-1,-1,2065	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,-1,-1,-1870	,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,-1,-1,1828	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,0,1,-1794	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,0,3,-1749	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,-1,1,-1565	,EP,D,M,MS,F,DB);
    DB = ADDdb(1,0,0,1,-1491	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,1,1,-1475	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,1,-1,-1410	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,0,-1,-1344	,EP,D,M,MS,F,DB);
    DB = ADDdb(1,0,0,-1,-1335	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,3,1,1107	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,0,-1,1021	,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,-1,1,833	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,1,-3,777	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,-2,1,671	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,0,-3,607	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,2,-1,596	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,1,-1,491	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,-2,1,-451	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,3,-1,439	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,2,1,422	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,0,-3,-1,421	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,1,-1,1,-366	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,1,0,1,-351	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,0,1,331	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,1,1,315	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-2,0,-1,302	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,0,1,3,-283	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,1,1,-1,-229	,EP,D,M,MS,F,DB);
    DB = ADDdb(1,1,0,-1,223	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(1,1,0,1,223	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,-2,-1,-220	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,1,-1,-1,-220	,EP,D,M,MS,F,DB);
    DB = ADDdb(1,0,1,1,-185	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(2,-1,-2,-1,181	,EP,D,M,MS,F,DB);
    DB = ADDdb(0,1,2,1,-177	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,-2,-1,176	,EP,D,M,MS,F,DB);
    DB = ADDdb(4,-1,-1,-1,166	,EP,D,M,MS,F,DB);
    DB = ADDdb(1,0,1,-1,-164	,EP,D,M,MS,F,DB);
    DB = ADDdb(4,0,1,-1,132	    ,EP,D,M,MS,F,DB);
    DB = ADDdb(1,0,-1,-1,-119	,EP,D,M,MS,F,DB);
    DB = ADDdb(4,-1,0,-1,115	,EP,D,M,MS,F,DB);
    DB = ADDdb(2,2,0,1,107	    ,EP,D,M,MS,F,DB);   

  % Zusatzterm für Laenge & Breite
  DLR(1) = DLR(1) + 3958 * sind(A1) + ...
           1962 * sind(L0-F) + 318 *sind(A2);
  DB = DB - 2235 * sind(L0) + ...
           382 * sind(A3) + 175 *sind(A1 - F) + ...
           175 * sind(A1 + F) + 127 * sind(L0 - MS) - ...
           115 * sind(L0 + MS);

  la  = L0 + DLR(1)/1e6;
  del = 385000.56 + DLR(2)/1000;    
  ba  = DB/1e6;
  
% Abstand in verschiedenen Masseinheiten
  if strcmpi(rMass,'km') 
       fac =1;  % r in km
  else
      if strcmpi(rMass,'RE') 
         fac = 1/6378.137; % r in Erdradien
      else
         fac =1/149597870.700;  % r in AE Einheiten
      end
  end
  MoonPos.ekl = [del*fac; deg2rad(la); deg2rad(ba)]; % geo-ekliptikal
  MoonPos.xyz = CalcXYZfromAngles(MoonPos.ekl); % xyz geo ekl 
  MoonPos.geo = mtimes(R_x(-epsE),MoonPos.xyz); % xyz in Equ-System  
  MoonPos.equ = CalcAnglesfromXYZ(MoonPos.geo); % RA, DEC 
end


function  DLR =ADDdldr(P, Q, R, S, CoefL, CoefR, EP, D, M, MS, F, DLR)
    DLRadd(1) =  CoefL*sind(P*D+Q*M+R*MS+S*F);
    DLRadd(2) =  CoefR*cosd(P*D+Q*M+R*MS+S*F);
    Qabs = abs(Q);
    if Qabs > 0 
       DLRadd(1) = DLRadd(1)*EP;
       DLRadd(2) = DLRadd(2)*EP;
    end
    if Qabs > 1 
       DLRadd(1) = DLRadd(1)*EP;
       DLRadd(2) = DLRadd(2)*EP;
    end
    DLR(1) = DLR(1) + DLRadd(1);
    DLR(2) = DLR(2) + DLRadd(2);
end

function  DB =ADDdb(P, Q, R, S, CoefB, EP, D, M, MS, F, DB)
    DBadd =  CoefB*sind(P*D+Q*M+R*MS+S*F);
%     Qabs = abs(Q);
%     if Qabs > 0 
%        DBadd = DBadd*EP;
%     end
%     if Qabs > 1 
%        DBadd = DBadd*EP;
%     end
    DB = DB + DBadd;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------