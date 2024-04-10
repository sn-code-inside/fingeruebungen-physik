% -------------------------------------------------------------------------
% LambertICBM.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Parameter des Abfangens einer ICBM.
% (Intercontinental Ballistic Missile)
% über die Lösung des Lambert-Problems 
%
% -------------------------------------------------------------------------

% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
GE  = 398600.4415;      % G*ME in km^3/s^2 fÃ¼r Erde
RE  = 6378;

%% Parameter Rechnung
% TOF
TOF = 1800; %in s
%Anfangsdaten Abfangrakete
vec_r1 =[6000, 3500, 0]';     % in km 
vec_v1 = [-2.5, 6.5, 2.5]';   % in km/s
r1     = vecnorm(vec_r1);     % in km 
v1     = vecnorm(vec_v1);     % in km/s

%Anfangsdaten ICBM
vec_R1 = [12250.0, 10250.0, 2000]'; % in km 
vec_V1 = [-3.5, 1.0, 0]';           % in km/s
R1     = vecnorm(vec_R1);           % in km 
V1     = vecnorm(vec_V1);           % in km/s


%% Bestimmung Bahn ICBM---------------------------------------------------% 
% 
% Wir benutzen die Rechnungen aus BahnparBestimmung01.m
BaPa = BahnPara(vec_R1,vec_V1,R1, V1,GE, TOF);
%Lösung Kepler-Gleichung um E2 und ups2 zu bestimmen, damit vec_R2 = vec_r2
%Alternative: Bei lleinen exz Nutzung der Mittelpunktsgleichung
%Alternative: Nutzung des Methodik Gauss-Vektoren aus der Himmelsmechanik

SatData = SatPQR(TOF, BaPa, GE);
vec_R2 = -SatData.xyz;          % in km 
R2     = vecnorm(vec_R2);       % in km 

% Trajektorie
tspan  = linspace(0,2*TOF,2*1801);
tspanB = linspace(0,-6*TOF,2*1801);
for  k=1:length(tspan)
   SatData = SatPQR(tspan(k), BaPa, GE);
   PosTarget(k,:) = -SatData.xyz; 
   if vecnorm(PosTarget(k,:)) < RE 
       PosTarget(k,:) = NaN;
   end
   SatData = SatPQR(tspanB(k), BaPa, GE);
   PosTargetB(k,:) = -SatData.xyz;    
   if vecnorm(PosTargetB(k,:)) < RE 
       PosTargetB(k,:) = NaN;
   end
end


%% Bestimmung Bahn Abfänger-----------------------------------------------% 
% 
% Wir benutzen die Rechnungen aus BahnparBestimmung01.m
BaPa = BahnPara(vec_r1,vec_v1,r1,v1,GE,TOF);
SatData = SatPQR(TOF, BaPa, GE);
vec_r2 = -SatData.xyz;                    % in km 
r2     = vecnorm(vec_r2);                 % in km 

% Trajektorie
tspan  = linspace(0,2*TOF,3601);
tspanB = linspace(0,-4*TOF,3601);
for  k=1:length(tspan)
   SatData = SatPQR(tspan(k), BaPa, GE);
   PosAbf(k,:) = -SatData.xyz;  
   if vecnorm(PosAbf(k,:)) < RE 
       PosAbf(k,:) = NaN;
   end
   SatData = SatPQR(tspanB(k), BaPa, GE);
   PosAbfB(k,:) = -SatData.xyz; 
   if vecnorm(PosAbfB(k,:)) < RE 
       PosAbfB(k,:) = NaN;
   end
end

%% Lambert-Problem

c     = vecnorm(vec_r1-vec_R2);
s     = (c+r1+R2)/2;

% Minimale Flugzeit (Parabel)
tmin     = sqrt(2)/3 * sqrt(s^3/GE) *(1-((s-c)/s)^3/2);
tmin_min = tmin/60;

% Maximale Flugzeit
amin     = s/2;
alpham   = 2*asin(sqrt(s/2/amin));
alphamd  = rad2deg(alpham);
betam    = 2*asin(sqrt((s-c)/2/amin));
betamd   = rad2deg(betam);
tmax     = sqrt(amin^3/GE)*(alpham-betam-(sin(alpham)-sin(betam)));
tmax_min = tmax/60;

% Iteration/Bisektion
[af, alphaf, betaf] = BisektionLambert(c, s, TOF, GE);

% Berechnung der notwendigen Geschwindigkeit des Abfängers am Punkt P1
A = sqrt(GE/4/af)*cot(alphaf/2);
B = sqrt(GE/4/af)*cot(betaf/2);
vec_ec = (vec_R2-vec_r1)/c;
vec_e1 = vec_r1/r1;
vec_v_tr1 = (B+A)*vec_ec + (B-A)*vec_e1;
v_tr1     = vecnorm(vec_v_tr1);

% Berechnung des Antriebsbedarfs
vec_Deltav = vec_v_tr1 - vec_v1;
Deltav     = vecnorm(vec_Deltav);

% Berechnung der Trajektorie des Abfängers nach Delta v
BaPa     = BahnPara(vec_r1,vec_v_tr1,r1,v_tr1, GE,TOF);
SatData  = SatPQR(TOF, BaPa, GE);
vec_r2tr = -SatData.xyz;                    % in km 
r2tr     = vecnorm(vec_r2tr);               % in km 

% Neue Trajektorie
tspan  = linspace(0,TOF,1801);
for  k=1:length(tspan)
   SatData       = SatPQR(tspan(k), BaPa, GE);
   PosAbftr(k,:) = -SatData.xyz;  
end


%% Graphik

figure(1)
AxEnd =20000;
hold on
view(3)
[x,y,z] = sphere;
surf(RE*x,RE*y,RE*z,'FaceAlpha',0.85);
axis([-AxEnd AxEnd -AxEnd AxEnd -AxEnd AxEnd]); 
xlabel('x in km'); ylabel('y in km'); zlabel('z in km')
plot3(PosTarget(:,1),PosTarget(:,2),PosTarget(:,3),...
      'LineWidth',2,'Color',Colors(4,:));
plot3(PosTargetB(:,1),PosTargetB(:,2),PosTargetB(:,3),...
      'LineWidth',2,'Color',Colors(4,:),'LineStyle',Style(4));
plot3(PosTarget(1,1),PosTarget(1,2),PosTarget(1,3),'d',...
      'MarkerSize',10,'MarkerFaceColor',Colors(4,:));
plot3(vec_R2(1),vec_R2(2),vec_R2(3),'s',...
      'MarkerSize',10,'MarkerFaceColor',Colors(4,:));
plot3(PosAbf(:,1),PosAbf(:,2),PosAbf(:,3),...
      'LineWidth',2,'Color',Colors(9,:));
plot3(PosAbfB(:,1),PosAbfB(:,2),PosAbfB(:,3),...
      'LineWidth',2,'Color',Colors(9,:),'LineStyle',Style(4));
plot3(PosAbf(1,1),PosAbf(1,2),PosAbf(1,3),'d',...
      'MarkerSize',10,'MarkerFaceColor',Colors(9,:));
plot3(vec_r2(1),vec_r2(2),vec_r2(3),'s',...
      'MarkerSize',10,'MarkerFaceColor',Colors(9,:));
plot3(PosAbftr(:,1),PosAbftr(:,2),PosAbftr(:,3),...
      'LineWidth',2,'Color',Colors(1,:));
axis equal
grid on

%-------------------------------------------------------------------------%
% Ende Programm
% ------------------------------------------------------------------------%

% Bahnparameter
function BaPa = BahnPara(vec_R1,vec_V1,R1, V1,GE, t)
    hv      = cross(vec_R1,vec_V1);                 % Drehimpuls
    h       = sqrt(dot(hv,hv));
    Wv      = hv/h;
    inkl    = atand(sqrt(hv(1)^2+hv(2)^2)/hv(3));   % Inklination
    Omega   = atan2d(Wv(1),-Wv(2));                 % Knoten
    p       = h^2/GE;                               % semi-latus rectum
    a       = 1/(2/R1-V1^2/GE);                     % groÃe Halbachse
    nBew    = sqrt(GE/a^3);                         % mittl. Bewegung, 3. KG
    exz     = sqrt(1-p/a);                          % ExzentrizitÃ¤t
    E1      = atan2(dot(vec_R1,vec_V1)/(a^2*nBew),(1-R1/a)); % Exzentr. Anomalie
    M1      = E1-exz*sin(E1);                       % Mittl. Anomalie   
    M1d     = rad2deg(M1);                          % Mittl. Anomalie
    u       = atand(vec_R1(3)/(-vec_R1(1)*Wv(2)+vec_R1(2)*Wv(1)));  % LÃ¤nge
    ups1    = atand(sqrt(1-exz^2)*sin(E1)/(cos(E1)-exz));   % wahre Anomalie
    omega   = u-ups1;                               % Argument Perihel
    M2      = deg2rad(M1d) + nBew*t;
    M2d     = rad2deg(M2);
    BaPa    = [a exz inkl M1d Omega omega];
end


% Iteration/Bisektion
function [a, alpha, beta] = BisektionLambert(c, s, TOF, GE)
    amax    = 2*s;
    amin    = s/2;
    n       = 1;
    aiter   = zeros(1,100);
    tau_n   = zeros(1,100);
    aiter(n)= (amin+amax)/2;
    while abs(1-tau_n(n)/TOF) > 0.001 
        alpha_n(n)  = 2*asin(sqrt(s/2/aiter(n)));
        beta_n(n)   = 2*asin(sqrt((s-c)/2/aiter(n)));
        tau_n(n+1)   = sqrt(aiter(n)^3/GE)*...
                   (alpha_n(n)-beta_n(n)-(sin(alpha_n(n))-sin(beta_n(n))));
        if tau_n(n+1) > TOF 
            amin = aiter(n);
        else
            amax = aiter(n);
        end
        n = n+1;
        aiter(n)= (amin+amax)/2;
    end
    nend  = n-1;
    a     = aiter(nend);
    alpha = alpha_n(nend);
    beta  = beta_n(nend);
end


function SatData= SatPQR(T, BaPa, GE) 
  % Zeit T in s, alle Winkel in deg
  aPn = BaPa(1);
  ePn = BaPa(2);
  iPn = BaPa(3);
  MPn = BaPa(4);  %Mittlere Anomalie zum Startzeitpunkt in deg
  OmegaPn = BaPa(5);
  omegaPn = BaPa(6);
  nP  = sqrt(GE/aPn^3);  %mittlere Bewegung in 1/s
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
  Vel     = sqrt(GE*aPn)./dis;
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
  SatData  = Sat;
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
