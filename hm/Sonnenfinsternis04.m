% -------------------------------------------------------------------------
% Sonnenfinsternis02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial basiert in Teilen 
% auf C++ Programmstrukturen und der Beschreibung/Erlaeuterung in  
%
% "Astronomie mit dem Personal Computer"
%
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Startet von berechneten ungefaehren Neumondzeiten und
% bestimmt, wann und wo um die Neumondzeit eine Sonnenfinsternis 
% stattfinden kann.
% Fuer die SoFi wird die Zentallinie berechnet.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
SunData = PerturbImport('SonnePos.csv');

% Bestimmung der Zentrallinien 
d_0(1,:)   = "08-Apr-2024 18:00:00"; %(UT)
d_0(2,:)   = "12-Aug-2026 17:00:00"; %(UT)
d_0(3,:)   = "02-Aug-2027 08:00:00"; %(UT)
d_0(4,:)   = "26-Jan-2028 13:00:00"; %(UT)
d_1(1,:)   = "08.04.2024"; %(UT)
d_1(2,:)   = "12.08.2026"; %(UT)
d_1(3,:)   = "02.08.2027"; %(UT)
d_1(4,:)   = "26.01.2028"; %(UT)

jend =4;

% Schleife ueber 4 SoFi-Zeitpunkte
for j=1:jend
    dt_0  = datetime(d_0(j,:),'InputFormat','dd-MMM-yyyy HH:mm:ss');
    T_0   = juliandate(dt_0)+0.000001; % UT
    DelTSec = ETminusUT(T_0);
    eps0=deg2rad(EpsErde(T_0));
    T_F   = T_0 + DelTSec/86400; % ET
    t_F   = Jd2JJht(T_F);  % Zeit in Julianischen Jahrhunderten seit JD2000
    delt = 2/24/60/36525; % Schrittweite in min 
    t_Begin = t_F-60*delt;
    t_End  = t_F+120*delt;
    StepMarkers = 30; % Marker jede Stunde fuer Graphik
    t = t_Begin; % ET in Jht
    k = 1;
    l = 1;
    fprintf(' \n');
    while t < t_End
        T = JJht2Jd(t);
        T_UT= T - (DelTSec)/86400;
        Zeit =   datetime(T_UT,'ConvertFrom','juliandate');
        ZL = Central (t, DelTSec, SunData);
        outstr = string(Zeit,'dd-MMM-yyyy HH:mm')+" UT";
        if ZL.Phase == 'total' || ZL.Phase == 'ringfoermig'
            lat(j,k) = ZL.geo(3);
            lon(j,k) = ZL.geo(2);
        else 
            lat(j,k) = NaN;
            lon(j,k) = NaN;
        end     
        fprintf(' %s   %s°   %s  %s   %s \n',...
                 outstr, StrDMS(ZL.geo(3)),StrDMS(ZL.geo(2)),...
                 StrHMS(ZL.T_KS/60), ZL.Phase);
        if mod(k,StepMarkers)- 1 ==0 %&& ZL.Phase == 'total'
               outstr1 = string(Zeit,'HH:mm')+" UT";
               HourLabel(j).Time(l)= k; 
               HourLabel(j).Label(l)= outstr1; 
               l = l+1;
        end
        t = t + delt;
        k = k+1;      
    end
end

%_________________________________________________________________________
%%

% Graphische Ausgabe
figure();
geoplot(lat(1,:),lon(1,:),'-+','MarkerIndices',HourLabel(1).Time(:),...
        'MarkerSize',8,'Color', Colors(4,:),'LineWidth', 2);
hold on
if jend > 1
    geoplot(lat(2,:),lon(2,:),'-+','MarkerIndices',HourLabel(2).Time(:),...
     'MarkerSize',8,'LineWidth', 2,'Color', Colors(10,:));
    geoplot(lat(3,:),lon(3,:),'-+','MarkerIndices',HourLabel(3).Time(:),...
     'MarkerSize',8,'Color', Colors(4,:),'LineWidth', 2,'LineStyle',':');
    geoplot(lat(4,:),lon(4,:),'-+','MarkerIndices',HourLabel(4).Time(:),...
     'MarkerSize',8,'Color', Colors(10,:),'LineWidth', 2,'LineStyle',':');
end
geobasemap('colorterrain');
for j =1:jend
    for l=1:length(HourLabel(j).Time)
        p=text(lat(j,HourLabel(j).Time(l)),lon(j,HourLabel(j).Time(l)),...
               HourLabel(j).Label(l),'VerticalAlignment','bottom',...
               'HorizontalAlignment','left');
        p.FontSize =14;
    end
end
if jend > 1
    legend (d_1,'location','southeast');
    legend box off;
end
geolimits([-60 60],[-180 +180])
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%__________________________________________________________________________
%%
%  Funktionen
%--------------------------------------------------------------------------                                                                           
% Central                                                                                                                                               
% Gibt die Zentrallinie in geographischen Koordinaten aus inkl. Phase der 
% Finsternis
function ZL = Central (t, DelTSec, SunData)
   % Werte brauchen wir zur Berechnung der Schattengeschwindigkeit
   dt = 0.1;                % kleines Zeitintervall in Minuten
   MPC= 52596000;           % Minuten pro Julian. Jahrhundert
   Omega = 4.3755e-3;       % Winkelgeschwindigkeit der Erde in rad/min
   tau = 8.32/1440/36525;   % Lichtlaufzeit Sonne-Erde;
   t_UT = t - DelTSec/86400/36525;
   epsE = deg2rad(EpsErde(JJht2Jd(t)));
   
   % Koordinatenberechnung
   MoonPos = MondExakt(t,epsE,'RE');
   SunPos  = SonneExakt(t-tau,SunData,epsE,'RE'); % Beachte tau ! 
   % Berechnung der Schattenachse der SoFi
   SA = SAchse(MoonPos, SunPos); 
   if SA.Phase == "total" || SA.Phase == "ringfoermig"
        MJD_UT   = t_UT*36525+51544.5;
        Dec = SA.equ(3);
        RA  = rad2deg(SA.equ(2));
        phi      = rad2deg(Dec) + 0.1924*sin(2*Dec); %in rad
        lambda   = wrapTo180(-15*LMST(MJD_UT,0.0)+RA); 
        MoonPos2 = MondExakt(t+dt/MPC,epsE,'RE');
        SunPos2  = SonneExakt(t+dt/MPC-tau,SunData,epsE,'RE'); 
        SA2 = SAchse(MoonPos2, SunPos2);
        if SA2.Phase == "ringfoermig" 
            MoonPos2 = MondExakt(t-dt/MPC,epsE,'RE');
            SunPos2  = SonneExakt(t-dt/MPC-tau,SunData,epsE,'RE'); 
            SA2 = SAchse(MoonPos2, SunPos2);
        end
        W = dt*Omega;
        DX = SA2.xyz(1)-SA.xyz(1)+W*SA.xyz(2);
        DY = SA2.xyz(2)-SA.xyz(2)-W*SA.xyz(1);
        DZ = SA2.xyz(3)-SA.xyz(3);
        VecD = [DX;DY;DZ];
        D    = sqrt(norm(VecD)^2-sum(VecD .*SA.E_XS)*sum(VecD .*SA.E_XS));
        T_KS = dt*abs(SA.DKS)/D;       
   else
        lambda = 0;
        phi    = 0; 
        t_KS   = 0; 
        T_KS   = 0;
   end
   ZL.geo   = [1; lambda; phi];
   ZL.Phase = SA.Phase;
   ZL.T_KS  = T_KS;
end
% -------------------------------------------------------------------------                                                                           
% Ende Funktion Central
% -------------------------------------------------------------------------                                                                           



%%
% -------------------------------------------------------------------------                                                                           
% Schattenachse SAchse                                                                                                                                               
% Berechnung der Schattenachse der SoFi
function SA = SAchse(MoonKoord, SunKoord)
   fac = 0.996633; % Abplattung Erde
   DM  = 0.5450;   % Monddurchmesser in Erdradien
   DS  = 218.25;   % Sonnendurchmesser in Erdradien
   X_M = MoonKoord.geo; % kartes. aequatoriale Koord. v. Mond und Sonne
   X_S = SunKoord.geo;
   X_M(3) = X_M(3)/fac; % Streckung der z-Koordinate wg Abplattung
   X_S(3) = X_S(3)/fac; 
   X_MS   = X_M-X_S;
   R_MS   = norm(X_MS);
   E_XS   = X_MS/R_MS;          % Richtungsvektor Sonne-Mond
   p0     = -sum(X_M .* E_XS);  % Entfernung Mond-Hauptebene
   Dqua   = p0*p0 + 1.0-norm(X_M)^2;
   r0     = sqrt(1-Dqua);       % Entfernung Erdmitte-Schattenachse
   DKS    = (DS-DM)*p0/R_MS-DM; % Durchmesser Kernschatten auf Hauptebene
   DHS    = (DS+DM)*p0/R_MS+DM; % Durchmesser Halbschatten auf Hauptebene
   XB = [NaN;NaN;NaN];
   if r0 < 1.0 %Schattenachse trifft auf Erde
      p   = p0 -sqrt(Dqua);
      DKS = (DS-DM)*(p/R_MS)-DM;% Durchmesser Kernschatten auf Erdoberflaeche
      XB  = X_M + p*E_XS;       % kartesisches Koordinaten des Schattenkegels
      XB(3) = XB(3)*fac;        % Rueckrechnung der Abplattung
      if DKS > 0 
          Phase = "ringfoermig";
      else
          Phase = "total";
      end
   else
      if r0 < (1.0+0.5*abs(DKS))
         if DKS > 0 
             Phase = "ringfoermig, nicht zentral";
         else 
             Phase = "total, nicht zentral";
         end
      else
          if r0 < (1.0+0.5*abs(DHS))
             Phase = "partiell";
          else
             Phase = "keine Finsternis";
          end
      end 
   end
   SA.Phase = Phase;
   SA.xyz   = XB;
   SA.equ   = CalcAnglesfromXYZ(SA.xyz);
   SA.E_XS  = E_XS;
   SA.DKS   = DKS;
end
% -------------------------------------------------------------------------                                                                           
% Ende Funktionen
% -------------------------------------------------------------------------                                                                           

