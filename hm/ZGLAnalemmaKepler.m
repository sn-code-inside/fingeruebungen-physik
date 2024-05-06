% -------------------------------------------------------------------------
% ZGLAnalemmaKepler.m
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
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Programm berechnet die planetenzentrischen Positionen der
% Sonne von den verschiedenen Planeten aus 
% auf Basis der Koordinaten aus der Keplergleichung und der 
% Stoerungstheorie (Venus bis Jupiter).
% Daraus wird die ZGL und das Analemma berechnet.
% Die Nullmeridianeichung wird o.B.d.A. auf den 1.1.2000 gelegt.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
AE      = 149597870;
c_light   = 1/0.00578;   % AE/d

dt1 = datetime('2000-01-01 12:00:00');  % UTC
T1  = juliandate(dt1);                  % Julianisches Datum  UT
t1  = Jd2JJht(T1);
T1E = T1 + ETminusUT(T1)/86400;         % Julianisches Datum  ET
t1E = Jd2JJht(T1E);
% Achsneigung Erde
epsErad = deg2rad(EpsErde(T1E));

% Mit oder ohne Light Time Correction
LTCorr = 0;  % mit = 1 ohne =  0

Aequi = 'J2000';  %Es spielt keine Rolle ob J2000 oder Datum, da die KOS
                  %fuer wahre und mittlere gleich sind und nur die
                  %Differenz eine Rolle spielt.
[BaPa,BaPadot]= OrbitParameter(T1,Aequi); 
EarthDay      = 86400; %Laenge eines Erdtages 

PlSelect1=3;  %Planetenauswahl Start Erde
PlSelect2=8;  %Planetenauswahl Ende  Neptun

% Eigentliche Berechnung
for PlNr = PlSelect1:PlSelect2
    Planets(PlNr)=PlanetPQR(T1E, BaPa, BaPadot, PlNr);
    %Verschieden Bahnparameter
    PiqP    = BaPa.PiqP(PlNr);
    Piqdot  = BaPadot.Piqdot(PlNr);
    n_anom  = BaPa.NP(PlNr);
    n_trop  = BaPadot.L0dot(PlNr);
    % Laenge eines Orbits/Jahres in Jul. Erdtagen
    NrJDpOrbit(PlNr) = 365.25*36000/n_trop;
    NrPoints         = 720;    %Stuetzstellen, Stuetzweite ist 0.5 Grad
    % Laenge eines Solartages (sol)
    SidTag           = BaPa.SidTagP(PlNr);     %Siderischer Tag
    SolTag           = SidTag/(1-SidTag/(NrJDpOrbit(PlNr)*86400));    
                                               %Solartag auf Planeten (sol)
    NrSolpOrb(PlNr)  = NrJDpOrbit(PlNr)*86400/SolTag;    
                                               %Solartage eines Jahres 
    
    % heliozentrischer Vektor zum Planeten in der Aequatorebene der Erde
    r       = mtimes(R_x(-epsErad),Planets(PlNr).xyz); 
    % Abstand Sonne-Planet
    Del     = norm(r);
    % Lichtlaufzeitkorrektur
    T1E_LT       = T1E - LTCorr*Del/c_light;
    t1E_LT       = Jd2JJht(T1E_LT);
    Planets(PlNr)= PlanetPQR(T1E_LT, BaPa, BaPadot, PlNr);

    % Berechnungsvektor ueber 720 Stuetzstellen
    TFE     = T1E   +NrJDpOrbit(PlNr);%+100;      % ET nach einem Umlauf
    TFE_LT  = T1E_LT+NrJDpOrbit(PlNr);%+100;      % ET_korr nach Umlauf
    tag     = linspace(T1E,TFE_LT,NrPoints+1); % Jul. ET, LT korrigiert
    TF      = TFE-ETminusUT(TFE)/86400;        % UT nach einem Umlauf
    tagUT   = linspace(T1,TF,NrPoints+1);      % Jul. UT   
    
    % Schleife fuer Umlauf
    for m=1:NrPoints+1
        dtag(PlNr,m) = datetime(tag(m),'convertfrom','juliandate');
        dtUT(PlNr,m) = datetime(tagUT(m),'convertfrom','juliandate');
        %heliozentrische Planetenposition
        Planets(PlNr)=PlanetPQR(tag(m),BaPa,BaPadot,PlNr);
        r           = mtimes(R_x(-epsErad),Planets(PlNr).xyz);
        [E, Sense]  = Orient (PlNr, t1E_LT);
        %planetozentrische Sonnenposition
        [R_equ, f]  = Shape (PlNr);
        [lonS, latS, lat2S] = Rotation (-r,E,Sense,f);
        %planetozentrische RA der wahren Sonne
        RAS(PlNr,m)  = wrapTo360(-rad2deg(lonS));
        %planetozentrische RA der mittleren Sonne
        %reversal for Uranus
        if PlNr == 7
            RAMS(PlNr,m) = wrapTo360(PiqP - (Piqdot+ n_anom)*Jd2JJht(tag(m))); 
        else
            RAMS(PlNr,m) = wrapTo360(PiqP + (Piqdot+ n_anom)*Jd2JJht(tag(m))); 
        end;
        Dec(PlNr,m)  = rad2deg(latS);
        %Zeitgleichung
        ZGL(PlNr,m)  = -4*wrapTo360(RAS(PlNr,m)-RAMS(PlNr,m));
        
        %Umkehrung Vorzeichen Uranus
        if PlNr == 7
            ZGL(PlNr,m) = -ZGL(PlNr,m);
        end;
    end
end 

%Berechnungen Zeitgleichungsmittelwert fuer eine angenommene Kreisbahn
for PlNr = PlSelect1:PlSelect2
    Planets(PlNr)=PlanetPQR(T1E, BaPa, BaPadot, PlNr);
    %Verschieden Bahnparameter
    PiqP    = BaPa.PiqP(PlNr);
    n_anom  = BaPa.NP(PlNr);
    n_trop  = BaPadot.L0dot(PlNr);
    % heliozentrischer Vektor zum Planeten in der Aequatorebene der Erde
    r       = mtimes(R_x(-epsErad),Planets(PlNr).xyz); 
    % Abstand Sonne-Planet
    Del     = norm(r);
    % Lichtlaufzeitkorrektur
    T1E_LT       = T1E - LTCorr*Del/c_light;
    t1E_LT       = Jd2JJht(T1E_LT);
    %heliozentrische Planetenposition
    BaPa.eP(PlNr)  = 0;  %Ab hier wird Kreisbahn angenommen !
    Planets(PlNr)=PlanetPQR(T1E_LT,BaPa,BaPadot,PlNr);
    if PlNr > 3 
        r           = mtimes(R_x(-epsErad),Planets(PlNr).xyz);
    else 
        r           = Planets(PlNr).xyz; 
    end
    [E, Sense]  = Orient (PlNr,t1E_LT);
    %planetozentrische Sonnenposition
    [R_equ, f]  = Shape (PlNr);
    [lonS, latS, lat2S] = Rotation (-r,E,Sense,f);
    %planetozentrische RA der wahren Sonne
    RAS0(PlNr)  = wrapTo360(-rad2deg(lonS));
    %planetozentrische RA der mittleren Sonne
    RAMS0(PlNr) = wrapTo360(PiqP + n_anom*t1E); 
    %Zeitgleichungsmittelwert fuer eine angenomemne Kreisbahn
    ZGL0(PlNr)  = -4*wrapTo360(RAS0(PlNr)-RAMS0(PlNr));
    
    %durch Integration ermittelter Mittelwert !!! wichtig!!!
    integral(PlNr) = trapz(tag-T1E,ZGL(PlNr,:));
    offset(PlNr) = integral(PlNr)/(TFE_LT-T1E);
end

% Zeiteichung korrigiert um Mittelwert
for PlNr = PlSelect1:PlSelect2
   %ZGL(PlNr,:) = ZGL(PlNr,:) - ZGL0(PlNr);
   ZGL(PlNr,:) = ZGL(PlNr,:) - offset(PlNr);
end

%%
%%
% Graphische Ausgabe

% Aufbereitung zum Ausgeben
for PlNr=PlSelect1:PlSelect2
    for m = NrPoints+1:length(Dec(PlNr,:))
       ZGL(PlNr,m) = NaN;  
       Dec(PlNr,m) = NaN;
    end
end

% Bild Zeitgleichung
figure(); 
for PlNr =PlSelect1:PlSelect2 %Darstellung ab Erde bis Neptun
    subplot(2,(PlSelect2-2)/2,PlNr-2);
    plot(dtag(PlNr,:),ZGL(PlNr,:),...
        'Color', Colors(PlNr,:),'LineStyle','-','LineWidth',2);
    hold on;
    ylim([-30 +30]);
    if PlNr == 7
        ylim([-720 +720]);
    end
    if PlNr == 4 || PlNr == 6
        ylim([-60 +60]);
    end
    ylabel('Datum  UT')
    ylabel('ZGL in min');
    title(BaPa.Name(PlNr));
    grid minor;
    set(gca,'Fontsize',16);
end

% Bild Analemma
figure();
for PlNr =PlSelect1:PlSelect2
    subplot(2,(PlSelect2-2)/2,PlNr-2);
    if PlNr == 7
       plot(ZGL(PlNr,:)/60,Dec(PlNr,:),...
        'Color', Colors(PlNr,:),'LineStyle','-','LineWidth',2);
     else 
        plot(ZGL(PlNr,:),Dec(PlNr,:),...
        'Color', Colors(PlNr,:),'LineStyle','-','LineWidth',2);
    end
    grid on,
    grid minor;
    ylim([-30 30]);
    xlim([-60 60]);
    title(BaPa.Name(PlNr));
    ylabel('Deklination')
    xlabel('ZGL in min');
    if PlNr == 7
        ylim([-100 100]);
        xlim([-12 +12]);
        xlabel('ZGL in h');
    end
    set(gca,'Fontsize',16);
end

% Druckausgabe
fprintf('\n %s Keplerlvsung \n', string(dt1) );
fprintf('Planet \n \n');
for PlNr =PlSelect1:PlSelect2
    %Laenge eines Orbits/Jahres in Julianischen Tagen 
    n_trop= BaPadot.L0dot(PlNr);
    NrJDpOrbit = 365.25*36000/n_trop;           %Erdtage eines Jahres 

    %Laenge eines Sonnentages (sol) und Laenge eines Orbits/Jahres in sols
    SidTag        = BaPa.SidTagP(PlNr);         %Siderischer Tag
    SolTag        = SidTag/(1-SidTag/(NrJDpOrbit*86400)); 
                                                %Solartag auf Planeten(sol)
    NrSolDpOrbit  = NrJDpOrbit*86400/SolTag;    %Solartage eines Jahres 
    SolSidRatio   = SolTag/SidTag;              %Verhaeltnis Sol/Sid Tag Mars
    SolJDRatio    = SolTag/EarthDay;            %Verhaeltnis Sol/Erdtag
    SidJDRatio    = SidTag/EarthDay;            %Verhaeltnis Sid/Erdtag

    fprintf('%s: \t Sol: %s  \t \t Sid.Tag: %s  \t   Trop.Jahr: %10.4f in Erdtagen \t \t  Trop.Jahr: %10.4f in sols \n',...
               string(BaPa.Name(PlNr)), StrHMS(SolTag*24/86400),... 
               StrHMS(SidTag*24/86400),... 
               NrJDpOrbit,...
               NrJDpOrbit*86400/SolTag); 
    fprintf('%s: \t Sol: %8.6f zu Erdtagen \t Sid.Tag: %8.6f zu Erdtagen \t  \n \n',...
               string(BaPa.Name(PlNr)), SolJDRatio,... 
               SidJDRatio);
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% 
% Funktionen
% --------------------------------------------------------------------------
% Shape: Liefert den Aequatorradius und die Abplattung der Planeten
% Eingabe:
%   PlNr  Identifiziert den Planeten 
% Ausgabe:
%   R_equ     Aequatorradius in [km]
%   f         Geometrische Abplattung
%--------------------------------------------------------------------------

function [R_equ, f] = Shape (PlNr)
 switch (PlNr) 
    case 10    
        R_equ = 696000.0; f = 0.0;       
    case 1     
        R_equ =  2439.0; f = 0.0;       
    case 2     
        R_equ =  6051.0; f = 0.0;        
    case 3
        R_equ =  6378.14; f = 0.00335281; 
    case 4
        R_equ =  3393.4;  f = 0.0051865;  
    case 5
        R_equ = 71398.0;  f = 0.0648088;  
    case 6
        R_equ = 60000.0;  f = 0.1076209;  
    case 7
        R_equ = 25400.0;  f = 0.030;      
    case 8
        R_equ = 24300.0;  f = 0.0259;     
    case 9
        R_equ =  1500.0;  f = 0.0;
 end
end

%--------------------------------------------------------------------------
% Orient: Berechnet die Orientierung des planetozentrischen Bezugssystems +
%         liefert die Rotationsrichtung
% Eingabe:
%   PlNr  Identifiziert den Planeten (s. APC_Planets.h)
%   t         Zeit in Julianischen Jahrhunderten seit J2000
% Ausgabe:
%   ET        Transformationsmatrix vom geozentrischen, auf den mittleren 
%             Aequator und das Aequinoktium J2000 bezogenen KOS
%             zum koerperfesten aequatorialen Bezugssystem
%   Sense     Rotationsrichtung (direkt oder retrograd)
%--------------------------------------------------------------------------
function [ET, Sense] = Orient (PlNr,t)
  % Rektaszension und Deklination der Rotationsachse bezogen auf
  % Aequator und Ekliptik zur Epoche J2000; Orientierung des
  % Nullmeridians
  d = t*36525;  % Tage seit J2000
  switch (PlNr)
    case 10
        RA  = 286.13; Dec =  63.87; W   =  84.182 +  14.1844000*d;  
    case 1
        RA  = 281.01  -   0.033*t;  Dec =  61.45  -   0.005*t;    
        W   = 329.68  +   6.1385025*d;  
    case 2
        RA  = 272.76; Dec =  67.16; W   = 160.20  -   1.4813688*d;  
    case 3
        RA  =   0.00  -   0.641*t;  Dec =  90.00  -   0.557*t;
        W   = 190.16  + 360.9856235*d;  
    case 4
        RA  = 317.681 -   0.108*t; Dec =  52.886 -   0.061*t;
        W   = 177.901 + 350.8919830*d;  
    case 5
        RA  = 268.05  -   0.009*t; Dec =  64.49  +   0.003*t;
        W =  67.10  + 877.900*d;  
    case 6
        RA  =  40.589 -   0.036*t; Dec =  83.537 -   0.004*t;  
        W   = 227.2037 + 844.3*d;
    case 7
        RA  = 257.311; Dec = -15.175;
        W   = 203.81  - 501.1600928*d;   
    case 8
        N   = deg2rad(357.85+52.316*t);    
        RA  = 299.36  + 0.70*sin(N); 
        Dec =  43.46  - 0.51*cos(N); 
        W   = 253.18  + 536.3128492*d - 0.48*sin(N);  
    case 9
        RA  = 313.02; Dec =   9.09;
        W   = 236.77  -  56.3623195*d;
  end
  W   = wrapTo360(W);
  RA  = deg2rad(RA); 
  Dec = deg2rad(Dec); 
  W   = deg2rad(W);  
  % Transformation vom mittleren Erdaequator und Fruehlingspunkt J2000
  % zum koerperfesten Bezugssystem von Aequator und Nullmeridian
  ET = PQR (W, pi/2-Dec, pi/2+RA);
  % Rotationsrichtung
  if ( PlNr==2 || PlNr==7 || PlNr==10 )
    Sense = +1;
  else
    Sense = -1;
  end
end 
 

%--------------------------------------------------------------------------
% Rotation: Berechnet die planetographischen Koordinaten
% Eingabe:
%   r        Planetozentrische Koordinaten der Erde 
%            bezogen auf Erdaequator 
%            und Aequinoktium J2000
%   E        Transformationsmatrix
%   Sense    Rotationsrichtung (direkt oder retrograd)
%   f        Geometrische Abplattung
% Ausgabe:
%   lon      Planetographische Laenge der Erde in [rad]
%   lat_g    Planetographische Breite der Erde in [rad]
%   lat_c    Planetozentrische Breite der Erde in [rad]
%--------------------------------------------------------------------------

function [lon, lat, lat2] =  Rotation (r, E, Sense, f)
  % Planetozentrische Breite und Laenge
  s     = mtimes(E,r);
  plako = CalcAnglesfromXYZ(s);
  lon   = Sense*plako(2); 
  lon   = Modulo(lon,2*pi);
  lat   = plako(3);
  lat2  = atan(tan(lat)/((1.0-f)*(1.0-f)));
end

