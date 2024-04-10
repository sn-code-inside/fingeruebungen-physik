% -------------------------------------------------------------------------
% MondExakt.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die geozentrischen Koordinaten des Mondes.
% Basis Brownsche Mondtheorie
% Alle Werte beziehen sich auf  Äquinoktium des Datums
%
%  t     : Zeit in Julianischen Jahrhunderten seit JD2000
%  epsE  : Ekliptikschiefe in Bogenmaß
%  rMass : Angabe für Masseinheit für r
%  Ausgabe:
%  MoonPos.ekl      Geozentrische Ekliptikale Korordinaten des Mondes
%  MoonPos.xyz      xyz Koordinaten des Mondes bezogen auf Erdmittelpunkt
%  MoonPos.geo      Geozentrische äquatoriale Korordinaten des Mondes
%  MoonPos.equ      Geozentrische äquatoriale xyz Koordinaten des Mondes 
% -------------------------------------------------------------------------

function MoonPos = MondExakt(t,epsE,rMass)  
pi2   = 2*pi;
Arcs  = 3600*180/pi;
t2    = t.*t; % Vereinfachung Zeit

CO = zeros(13,4);
SI = zeros(13,4);

% Stoerungen als Vektor auf Null bzw. Anfangswerte setzen
% ST(1)=DLAM=0; ST(2)=DS=0; ST(3)=GAM1C=0; ST(4)=SINPI=3422.7000;
ST  = [0 0 0 3422.7];
STN = 0.0;

% Langperiodische Stoerungen
S1 = Sine(0.19833+0.05611.*t);  S2 = Sine(0.27869+0.04508.*t);
S3 = Sine(0.16827-0.36903.*t);  S4 = Sine(0.34734-5.37261.*t);
S5 = Sine(0.10498-5.37899.*t);  S6 = Sine(0.42681-0.41855.*t);
S7 = Sine(0.14943-5.37511.*t); 

DL0 = 0.84.*S1+0.31*S2+14.27.*S3+ 7.26*S4+ 0.28*S5+0.24.*S6;
DL  = 2.94.*S1+0.31*S2+14.27.*S3+ 9.34*S4+ 1.12*S5+0.83.*S6;
DLS =-6.40.*S1                                    -1.89.*S6;
DF  = 0.21.*S1+0.31*S2+14.27.*S3-88.70.*S4-15.30*S5+0.24.*S6-1.86.*S7;
DD  = DL0-DLS;

DGAM   = -3332e-9 * Sine(0.59734-5.37261.*t) - ...
           539e-9 * Sine(0.35498-5.37899.*t) - ...
            64e-9 * Sine(0.39943-5.37511.*t);

% Mittlere Argumente der Mondbahn (inkl. langperiodische Korrekturen)
% L0 mittlere Laenge des Mondes
% l  mittlere Anomalie des Mondes  ls mittlere Anomalie der Sonne      
% F  mittlerer Knotenabstand       D  mittlere Elongation von der Sonne

L0 = pi2*Frac(0.60643382+1336.85522467*t-0.00000313*t2) + DL0/Arcs;
L  = pi2*Frac(0.37489701+1325.55240982*t+0.00002565*t2) + DL /Arcs;
LS = pi2*Frac(0.99312619+  99.99735956*t-0.00000044*t2) + DLS/Arcs;
F  = pi2*Frac(0.25909118+1342.22782980*t-0.00000892*t2) + DF /Arcs;
D  = pi2*Frac(0.82736186+1236.85308708*t-0.00000397*t2) + DD /Arcs;

for I=1:4
    os = 7;
    switch I
        case 1
            ARG = L;  MAX = 4; FAC = 1.000002208;
        case 2
            ARG = LS; MAX = 3; FAC = 0.997504612-0.002495388.*t;
        case 3
            ARG = F;  MAX = 4; FAC = 1.000002708+139.978*DGAM;
        otherwise
            ARG = D;  MAX = 6; FAC = 1.0;
    end
    % Index um 7 verschoben
    CO(os,I) = 1; CO(os+1,I) = cos(ARG)*FAC; CO(os-1,I)= +CO(os+1,I);
    SI(os,I) = 0; SI(os+1,I) = sin(ARG)*FAC; SI(os-1,I)= -SI(os+1,I);
    for J = 2:MAX
        [CO(os+J,I),SI(os+J,I)] = ADDTHE(CO(os+J-1,I),SI(os+J-1,I),...
                                  CO(os+1,I),SI(os+1,I));
        CO(os-J,I) = +CO(os+J,I);
        SI(os-J,I) = -SI(os+J,I);
    end
end

ST  = SOLAR1(ST,CO,SI);
ST  = SOLAR2(ST,CO,SI);
ST  = SOLAR3(ST,CO,SI);
STN = SOLARN(STN,CO,SI);
ST  = PLANETARY(t,ST);
% Laenge 
la  = 360*Frac((L0+ST(1)/Arcs)/pi2);
% Breite 
S   = F + ST(2)/Arcs;
FAC = 1.000002708+139.978*DGAM;
ba  = (FAC*(18518.511+1.189+ST(3))*sin(S)-6.24*sin(3*S)+STN)/3600.0;
% Distanz 
ST(4) = ST(4) * 0.999953253;
del   = Arcs/ST(4);
% Auswahl der Ausgabemoeglichkeit für Abstand
% Abstand in Erdradien
if strcmpi(rMass,'RE') 
   fac =1;  % r in Erdradien 
else
  if strcmpi(rMass,'km') 
     fac = 1*6378.137; % r in km
  else
     fac =149597870.700/6378.137;  % r in AE Einheiten
  end
end
MoonPos.ekl = [del*fac; deg2rad(la); deg2rad(ba)]; % geo-ekliptikal
MoonPos.xyz = CalcXYZfromAngles(MoonPos.ekl); % xyz geo ekl 
MoonPos.geo = mtimes(R_x(-epsE),MoonPos.xyz); % xyz in Equ-System  
MoonPos.equ = CalcAnglesfromXYZ(MoonPos.geo); % RA, DEC 
end


% Berechnet sin(x); in einheiten von pi2; x in rad
function SINE = Sine(x)
   SINE = sin(2*pi*Frac(x));
end


% Berechnet c=cos(a1+a2) und s=sin(a1+a2) aus den Additionstheoremen fuer 
% c1=cos(a1), s1=sin(a1), c2=cos(a2) und s2=sin(a2)   
function [C,S] = ADDTHE(c1,s1,c2,s2)
  C = c1*c2-s1*s2;
  S = s1*c2+c1*s2;
end

% Berechnet X=cos(p*arg1+q*arg2+r*arg3+s*arg4) und   
% Y=sin(p*arg1+q*arg2+r*arg3+s*arg4)                    

function [x,y]= TERM(P,Q,R,S,CO,SI)
    os = 7;    
    IZ = [P, Q, R, S];
    IZ = IZ+os;
    x = 1;
    y = 0;
    for k=1:4
        if ~(IZ(k) == os) 
             [x,y] = ADDTHE(x,y,CO(IZ(k),k),SI(IZ(k),k));
        end
    end
end

function ST = ADDSOL(COEFFL,COEFFS,COEFFG,COEFFP,P,Q,R,S,CO,SI,ST)
    [x,y]= TERM(P,Q,R,S,CO,SI);
    ST(1)= ST(1) + COEFFL*y; 
    ST(2)= ST(2) + COEFFS*y; 
    ST(3)= ST(3) + COEFFG*x; 
    ST(4)= ST(4) + COEFFP*x; 
end

function ST = SOLAR1(ST,CO,SI)
      ST = ADDSOL(   13.902,   14.06,-0.001,   0.2607,0, 0, 0, 4, CO,SI,ST);
      ST = ADDSOL(    0.403,   -4.01,+0.394,   0.0023,0, 0, 0, 3, CO,SI,ST);
      ST = ADDSOL( 2369.912, 2373.36,+0.601,  28.2333,0, 0, 0, 2, CO,SI,ST);
      ST = ADDSOL( -125.154, -112.79,-0.725,  -0.9781,0, 0, 0, 1, CO,SI,ST);
      ST = ADDSOL(    1.979,    6.98,-0.445,   0.0433,1, 0, 0, 4, CO,SI,ST);
      ST = ADDSOL(  191.953,  192.72,+0.029,   3.0861,1, 0, 0, 2, CO,SI,ST);
      ST = ADDSOL(   -8.466,  -13.51,+0.455,  -0.1093,1, 0, 0, 1, CO,SI,ST);
      ST = ADDSOL(22639.500,22609.07,+0.079, 186.5398,1, 0, 0, 0, CO,SI,ST);
      ST = ADDSOL(   18.609,    3.59,-0.094,   0.0118,1, 0, 0,-1, CO,SI,ST);
      ST = ADDSOL(-4586.465,-4578.13,-0.077,  34.3117,1, 0, 0,-2, CO,SI,ST);
      ST = ADDSOL(   +3.215,    5.44,+0.192,  -0.0386,1, 0, 0,-3, CO,SI,ST);
      ST = ADDSOL(  -38.428,  -38.64,+0.001,   0.6008,1, 0, 0,-4, CO,SI,ST);
      ST = ADDSOL(   -0.393,   -1.43,-0.092,   0.0086,1, 0, 0,-6, CO,SI,ST);
      ST = ADDSOL(   -0.289,   -1.59,+0.123,  -0.0053,0, 1, 0, 4, CO,SI,ST);
      ST = ADDSOL(  -24.420,  -25.10,+0.040,  -0.3000,0, 1, 0, 2, CO,SI,ST);
      ST = ADDSOL(   18.023,   17.93,+0.007,   0.1494,0, 1, 0, 1, CO,SI,ST);
      ST = ADDSOL( -668.146, -126.98,-1.302,  -0.3997,0, 1, 0, 0, CO,SI,ST);
      ST = ADDSOL(    0.560,    0.32,-0.001,  -0.0037,0, 1, 0,-1, CO,SI,ST);
      ST = ADDSOL( -165.145, -165.06,+0.054,   1.9178,0, 1, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -1.877,   -6.46,-0.416,   0.0339,0, 1, 0,-4, CO,SI,ST);
      ST = ADDSOL(    0.213,    1.02,-0.074,   0.0054,2, 0, 0, 4, CO,SI,ST);
      ST = ADDSOL(   14.387,   14.78,-0.017,   0.2833,2, 0, 0, 2, CO,SI,ST);
      ST = ADDSOL(   -0.586,   -1.20,+0.054,  -0.0100,2, 0, 0, 1, CO,SI,ST);
      ST = ADDSOL(  769.016,  767.96,+0.107,  10.1657,2, 0, 0, 0, CO,SI,ST);
      ST = ADDSOL(   +1.750,    2.01,-0.018,   0.0155,2, 0, 0,-1, CO,SI,ST);
      ST = ADDSOL( -211.656, -152.53,+5.679,  -0.3039,2, 0, 0,-2, CO,SI,ST);
      ST = ADDSOL(   +1.225,    0.91,-0.030,  -0.0088,2, 0, 0,-3, CO,SI,ST);
      ST = ADDSOL(  -30.773,  -34.07,-0.308,   0.3722,2, 0, 0,-4, CO,SI,ST);
      ST = ADDSOL(   -0.570,   -1.40,-0.074,   0.0109,2, 0, 0,-6, CO,SI,ST);
      ST = ADDSOL(   -2.921,  -11.75,+0.787,  -0.0484,1, 1, 0, 2, CO,SI,ST);
      ST = ADDSOL(   +1.267,    1.52,-0.022,   0.0164,1, 1, 0, 1, CO,SI,ST);
      ST = ADDSOL( -109.673, -115.18,+0.461,  -0.9490,1, 1, 0, 0, CO,SI,ST);
      ST = ADDSOL( -205.962, -182.36,+2.056,  +1.4437,1, 1, 0,-2, CO,SI,ST);
      ST = ADDSOL(    0.233,    0.36, 0.012,  -0.0025,1, 1, 0,-3, CO,SI,ST);
      ST = ADDSOL(   -4.391,   -9.66,-0.471,   0.0673,1, 1, 0,-4, CO,SI,ST);
end

function  ST = SOLAR2(ST,CO,SI)
      ST = ADDSOL(    0.283,    1.53,-0.111,  +0.0060,1,-1, 0,+4, CO,SI,ST);
      ST = ADDSOL(   14.577,   31.70,-1.540,  +0.2302,1,-1, 0, 2, CO,SI,ST);
      ST = ADDSOL(  147.687,  138.76,+0.679,  +1.1528,1,-1, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -1.089,    0.55,+0.021,   0.0   ,1,-1, 0,-1, CO,SI,ST);
      ST = ADDSOL(   28.475,   23.59,-0.443,  -0.2257,1,-1, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -0.276,   -0.38,-0.006,  -0.0036,1,-1, 0,-3, CO,SI,ST);
      ST = ADDSOL(    0.636,    2.27,+0.146,  -0.0102,1,-1, 0,-4, CO,SI,ST);
      ST = ADDSOL(   -0.189,   -1.68,+0.131,  -0.0028,0, 2, 0, 2, CO,SI,ST);
      ST = ADDSOL(   -7.486,   -0.66,-0.037,  -0.0086,0, 2, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -8.096,  -16.35,-0.740,   0.0918,0, 2, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -5.741,   -0.04, 0.0  ,  -0.0009,0, 0, 2, 2, CO,SI,ST);
      ST = ADDSOL(    0.255,    0.0 , 0.0  ,   0.0   ,0, 0, 2, 1, CO,SI,ST);
      ST = ADDSOL( -411.608,   -0.20, 0.0  ,  -0.0124,0, 0, 2, 0, CO,SI,ST);
      ST = ADDSOL(    0.584,    0.84, 0.0  ,  +0.0071,0, 0, 2,-1, CO,SI,ST);
      ST = ADDSOL(  -55.173,  -52.14, 0.0  ,  -0.1052,0, 0, 2,-2, CO,SI,ST);
      ST = ADDSOL(    0.254,    0.25, 0.0  ,  -0.0017,0, 0, 2,-3, CO,SI,ST);
      ST = ADDSOL(   +0.025,   -1.67, 0.0  ,  +0.0031,0, 0, 2,-4, CO,SI,ST);
      ST = ADDSOL(    1.060,    2.96,-0.166,   0.0243,3, 0, 0,+2, CO,SI,ST);
      ST = ADDSOL(   36.124,   50.64,-1.300,   0.6215,3, 0, 0, 0, CO,SI,ST);
      ST = ADDSOL(  -13.193,  -16.40,+0.258,  -0.1187,3, 0, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -1.187,   -0.74,+0.042,   0.0074,3, 0, 0,-4, CO,SI,ST);
      ST = ADDSOL(   -0.293,   -0.31,-0.002,   0.0046,3, 0, 0,-6, CO,SI,ST);
      ST = ADDSOL(   -0.290,   -1.45,+0.116,  -0.0051,2, 1, 0, 2, CO,SI,ST);
      ST = ADDSOL(   -7.649,  -10.56,+0.259,  -0.1038,2, 1, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -8.627,   -7.59,+0.078,  -0.0192,2, 1, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -2.740,   -2.54,+0.022,   0.0324,2, 1, 0,-4, CO,SI,ST);
      ST = ADDSOL(    1.181,    3.32,-0.212,   0.0213,2,-1, 0,+2, CO,SI,ST);
      ST = ADDSOL(    9.703,   11.67,-0.151,   0.1268,2,-1, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -0.352,   -0.37,+0.001,  -0.0028,2,-1, 0,-1, CO,SI,ST);
      ST = ADDSOL(   -2.494,   -1.17,-0.003,  -0.0017,2,-1, 0,-2, CO,SI,ST);
      ST = ADDSOL(    0.360,    0.20,-0.012,  -0.0043,2,-1, 0,-4, CO,SI,ST);
      ST = ADDSOL(   -1.167,   -1.25,+0.008,  -0.0106,1, 2, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -7.412,   -6.12,+0.117,   0.0484,1, 2, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -0.311,   -0.65,-0.032,   0.0044,1, 2, 0,-4, CO,SI,ST);
      ST = ADDSOL(   +0.757,    1.82,-0.105,   0.0112,1,-2, 0, 2, CO,SI,ST);
      ST = ADDSOL(   +2.580,    2.32,+0.027,   0.0196,1,-2, 0, 0, CO,SI,ST);
      ST = ADDSOL(   +2.533,    2.40,-0.014,  -0.0212,1,-2, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -0.344,   -0.57,-0.025,  +0.0036,0, 3, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -0.992,   -0.02, 0.0  ,   0.0   ,1, 0, 2, 2, CO,SI,ST);
      ST = ADDSOL(  -45.099,   -0.02, 0.0  ,  -0.0010,1, 0, 2, 0, CO,SI,ST);
      ST = ADDSOL(   -0.179,   -9.52, 0.0  ,  -0.0833,1, 0, 2,-2, CO,SI,ST);
      ST = ADDSOL(   -0.301,   -0.33, 0.0  ,   0.0014,1, 0, 2,-4, CO,SI,ST);
      ST = ADDSOL(   -6.382,   -3.37, 0.0  ,  -0.0481,1, 0,-2, 2, CO,SI,ST);
      ST = ADDSOL(   39.528,   85.13, 0.0  ,  -0.7136,1, 0,-2, 0, CO,SI,ST);
      ST = ADDSOL(    9.366,    0.71, 0.0  ,  -0.0112,1, 0,-2,-2, CO,SI,ST);
      ST = ADDSOL(    0.202,    0.02, 0.0  ,   0.0   ,1, 0,-2,-4, CO,SI,ST);
end  

function  ST = SOLAR3(ST, CO, SI)
      ST = ADDSOL(    0.415,    0.10, 0.0  ,  0.0013,0, 1, 2, 0, CO,SI,ST);
      ST = ADDSOL(   -2.152,   -2.26, 0.0  , -0.0066,0, 1, 2,-2, CO,SI,ST);
      ST = ADDSOL(   -1.440,   -1.30, 0.0  , +0.0014,0, 1,-2, 2, CO,SI,ST);
      ST = ADDSOL(    0.384,   -0.04, 0.0  ,  0.0   ,0, 1,-2,-2, CO,SI,ST);
      ST = ADDSOL(   +1.938,   +3.60,-0.145, +0.0401,4, 0, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -0.952,   -1.58,+0.052, -0.0130,4, 0, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -0.551,   -0.94,+0.032, -0.0097,3, 1, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -0.482,   -0.57,+0.005, -0.0045,3, 1, 0,-2, CO,SI,ST);
      ST = ADDSOL(    0.681,    0.96,-0.026,  0.0115,3,-1, 0, 0, CO,SI,ST);
      ST = ADDSOL(   -0.297,   -0.27, 0.002, -0.0009,2, 2, 0,-2, CO,SI,ST);
      ST = ADDSOL(    0.254,   +0.21,-0.003,  0.0   ,2,-2, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -0.250,   -0.22, 0.004,  0.0014,1, 3, 0,-2, CO,SI,ST);
      ST = ADDSOL(   -3.996,    0.0 , 0.0  , +0.0004,2, 0, 2, 0, CO,SI,ST);
      ST = ADDSOL(    0.557,   -0.75, 0.0  , -0.0090,2, 0, 2,-2, CO,SI,ST);
      ST = ADDSOL(   -0.459,   -0.38, 0.0  , -0.0053,2, 0,-2, 2, CO,SI,ST);
      ST = ADDSOL(   -1.298,    0.74, 0.0  , +0.0004,2, 0,-2, 0, CO,SI,ST);
      ST = ADDSOL(    0.538,    1.14, 0.0  , -0.0141,2, 0,-2,-2, CO,SI,ST);
      ST = ADDSOL(    0.263,    0.02, 0.0  ,  0.0   ,1, 1, 2, 0, CO,SI,ST);
      ST = ADDSOL(    0.426,   +0.07, 0.0  , -0.0006,1, 1,-2,-2, CO,SI,ST);
      ST = ADDSOL(   -0.304,   +0.03, 0.0  , +0.0003,1,-1, 2, 0, CO,SI,ST);
      ST = ADDSOL(   -0.372,   -0.19, 0.0  , -0.0027,1,-1,-2, 2, CO,SI,ST);
      ST = ADDSOL(   +0.418,    0.0 , 0.0  ,  0.0   ,0, 0, 4, 0, CO,SI,ST);
      ST = ADDSOL(   -0.330,   -0.04, 0.0  ,  0.0   ,3, 0, 2, 0, CO,SI,ST);
end

function  STN = ADDN(COEFFN,P,Q,R,S,STN,CO,SI)
      [x,y]= TERM(P,Q,R,S,CO,SI);
      STN = STN + COEFFN*y;
end

% Stoerungsanteil STN der ekliptikalen Breite           
function STN = SOLARN(STN, CO, SI)
%     STN = 0.0;
    STN = ADDN(-526.069, 0, 0,1,-2, STN,CO,SI);
    STN = ADDN(  -3.352, 0, 0,1,-4, STN,CO,SI);
    STN = ADDN( +44.297,+1, 0,1,-2, STN,CO,SI); 
    STN = ADDN(  -6.000,+1, 0,1,-4, STN,CO,SI);
    STN = ADDN( +20.599,-1, 0,1, 0, STN,CO,SI); 
    STN = ADDN( -30.598,-1, 0,1,-2, STN,CO,SI);
    STN = ADDN( -24.649,-2, 0,1, 0, STN,CO,SI); 
    STN = ADDN(  -2.000,-2, 0,1,-2, STN,CO,SI);
    STN = ADDN( -22.571, 0,+1,1,-2, STN,CO,SI); 
    STN = ADDN( +10.985, 0,-1,1,-2, STN,CO,SI);
end

% Stoerungen der ekliptikalen Laenge durch Venus und Jupiter           
function  ST = PLANETARY(t, ST)
   ST(1)  = ST(1) + ...
        0.82*Sine(0.7736  -62.5512.*t)+0.31*Sine(0.0466 -125.1025.*t) + ...
        0.35*Sine(0.5785  -25.1042.*t)+0.66*Sine(0.4591+1335.8075.*t) +...
        0.64*Sine(0.3130  -91.5680.*t)+1.14*Sine(0.1480+1331.2898.*t) +...
        0.21*Sine(0.5918+1056.5859.*t)+0.44*Sine(0.5784+1322.8595.*t) +...
        0.24*Sine(0.2275   -5.7374.*t)+0.28*Sine(0.2965   +2.6929.*t) +...
        0.33*Sine(0.3132   +6.3368.*t);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------