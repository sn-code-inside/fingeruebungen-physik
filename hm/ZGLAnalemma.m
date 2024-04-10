% -------------------------------------------------------------------------
% ZGLAnalemma.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Analemma von einem anderen Planeten
% 
% Das vorliegende Programmaterial basiert z.T. auf 
% Michael Allison und Megan McEwan
% Planetary and Space Science 48(2000) 215
%
% Programm berechnet die ZGL und das Analemma nach Kapitel 3.5.2.1.
% Die ZGL wird über
% ZGL = - C + (alpha_S - Laenge_S) berechnet.
% Der zweite Term ist die Reduktion auf den Aequator (siehe Beispiel 3.5.2)
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung Epoche
dt1 = datetime('2000-01-01 12:00:00');
T1  = juliandate(dt1);
t1  = Jd2JJht(T1);

epsE    = EpsErde(T1);  
epsErad = deg2rad(epsE);  

Aequi = 'J2000';  %Es pielt keine Rolle ob J2000 oder Datum, da die KOS
                  %fuer wahre und mittlere gleich sind und nur die
                  %Differenz eine Rolle spielt.
[BaPa,BaPadot]= OrbitParameter(T1,Aequi); 
EarthDay      = 86400; %Laenge eines Erdtages 


%%
%------------------------------------------------------------------------%%
% Beginn Rechnung

PlSelect1=3;  % Planetenauswahl Start Erde
PlSelect2=8;  % Planetenauswahl Ende  Neptun

for PlNr=PlSelect1:PlSelect2
    % Bahnparameter aus Datenfiles
    ex    = BaPa.eP(PlNr);          % Exzentrizitaet
    epsP  = BaPa.epsP(PlNr);        % Achsneigung
    M0P   = BaPa.M0P(PlNr);         % mittl. Anaomalie bei J2000
    n_anom= BaPa.NP(PlNr);          % anomal. Umlauf
    n_trop= BaPadot.L0dot(PlNr);    % trop. Umlauf
    M0dot = BaPadot.M0dot(PlNr);    % trop. Umlauf 
    L0P   = BaPa.L0P(PlNr);         % mittl. Laenge bei J2000
    PiqP  = BaPa.PiqP(PlNr);        % Laenge des Perihels bei J2000 
    Piqdot= BaPadot.Piqdot(PlNr);   % Praezession des Perihels 
    iP     = BaPa.iP(PlNr);         % orbital inclination mars
    OmegaP = BaPa.OmegaP(PlNr);     % longitude of node
    omega  = PiqP-OmegaP;           % angle of perihelion
    Omega90= OmegaP-90;
    if PlNr==7
        epsP = 82.23;  % Wert von Allison, 
     end

    if PlNr==1
        epsP = 8;  % Wert von Allison, 
     end

    % Berechnung der planetozentr. Laenge (!!!) des Perihels
    % bezogen auf Fruehlingspunkt des Planeten und der planetozentr.
    % "Rektaszension" der mittleren Sonne zur Epoche J2000

    [alpha0,delta0,Sense] = RotAxis(PlNr,Jd2JJht(T1));
    %convert to geocentric ecliptical coords
    axis_b = asind(-cosd(delta0)*sind(alpha0)*sind(epsE)+ ...
             sind(delta0)*cosd(epsE));
    axis_l = wrapTo360(atan2d(cosd(delta0)*sind(alpha0)*cosd(epsE)+ ... 
             sind(delta0)*sind(epsE),cosd(delta0)*cosd(alpha0))); 
    %rotation matrices
    Rz_O   = [cosd(OmegaP) -sind(OmegaP) 0; sind(OmegaP) cosd(OmegaP) 0; 0 0 1];
    Rz_mO  = [cosd(-OmegaP) -sind(-OmegaP) 0; sind(-OmegaP) cosd(-OmegaP) 0;0 0 1];
    Rz_O90 = [cosd(Omega90) -sind(Omega90) 0;...
             sind(Omega90) cosd(Omega90) 0;0 0 1];
    Rz_mO90= [cosd(-Omega90) -sind(-Omega90) 0;...
             sind(-Omega90) cosd(-Omega90) 0;0 0 1];
    Rz_o   = [cosd(omega) -sind(omega) 0; sind(omega) cosd(omega) 0;0 0 1];
    Rz_mo  = [cosd(-omega) -sind(-omega) 0; sind(-omega) cosd(-omega) 0;0 0 1];
    Rx_i   = [1 0 0;0 cosd(iP) -sind(iP); 0 sind(iP) cosd(iP)];
    Rx_mi  = [1 0 0;0 cosd(-iP) -sind(-iP); 0 sind(-iP) cosd(-iP)];
    
    %convert polar axis to cart. geocentric coords
    p      = [0,0,0];
    [p(1),p(2),p(3)] = sph2cart(axis_l/180*pi,axis_b/180*pi,1);
    
    %convert polar axis cart. geocentric coords to planeto-orbital cart. coords
    %can be achieved by two subsequent rotations of about z-axis with -omega 
    %and about x-axis with inclination -i
    p_bahn = Rx_mi*Rz_mO*p';
    
    %calculate angular momentum vector of orbit in orbital cart. coordinates
    L_bahn = [0,0,1];
    
    %calculate vernal point in orbit coordinates
    %vernal point is perpendicular to angular momentum and pependicular to
    %polar axis as it is the crossing point of equator and orbit
    vernal = cross(p_bahn',L_bahn);
  
    %convert vernal point orbital cart. coords to spherical coordinates
    [a,b,c]=cart2sph(vernal(1),vernal(2),vernal(3));
      
    %determine planeto-centric longitude of perihelion
    %by LSP = omega - a + 180° (180° because of reversal sun-planet coords
    LSP = wrapTo360(omega-rad2deg(a)+180);
    %determine right ascension of fictitious mean sun by
    % alpha_FMS = LSP + M0P
    alpha_FMS = LSP+M0P;

    % Faktoren fuer Mittelpunktsgleichung aus Exzentrizität
    fac(1) = (2*ex^1 - 1/4*ex^3   + 5/96*ex^5);
    fac(2) = (5/4*ex^2 - 11/24*ex^4 + 17/192*ex^6);
    fac(3) = (13/12*ex^3 - 43/63*ex^5);
    fac(4) = (103/96*ex^4 - 451/480*ex^6);
    fac(5) = (1097/960*ex^5);

    %Länge eines Orbits/Jahres in Julianischen Tagen 
    NrJDpOrbit = 365.25*36000/n_trop;           %Erdtage eines Jahres 

    %Länge eines Sonnentages (sol) und Länge eines Orbits/Jahres in sols
    SidTag        = BaPa.SidTagP(PlNr);         %Siderischer Tag
    SolTag        = SidTag/(1-SidTag/(NrJDpOrbit*86400));    %Solartag auf Planeten (sol)
    NrSolDpOrbit  = NrJDpOrbit*86400/SolTag;    %Solartage eines Jahres 
    SolSidRatio   = SolTag/SidTag;              %Verhältnis Sol/Sid Tag Mars
    SolJDRatio    = SolTag/EarthDay;            %Verhältnis Sol/Erdtag
    SidJDRatio    = SidTag/EarthDay;            %Verhältnis Sid/Erdtag
    % Julianische Startdatum für Berehcnungen in ET
    dt1 = datetime('2000-01-01 12:00:00');
    T1      = juliandate(dt1);

    %Datumsvektor für Zeitbereich
    NrDays(PlNr)= floor(NrSolDpOrbit)+1;                  %Sols für Berechnung
    Tv          = linspace(T1,T1+NrJDpOrbit+1,NrDays(PlNr)+1);  %Jul. Tage
    tv          = Jd2JJht(Tv);                            %Julian. Jahrhunderte 
    dtv(PlNr,1:length(Tv)) = datetime(Tv,'convertfrom','juliandate');

    % Schleife ueber die Julian. Tage eines Jahres
    for tag = 1:NrDays(PlNr)+1  %in Julianischen Tagen
        tz = tv(tag);
        M(tag)       = M0P + n_anom*tz;             %Mittl. Anomalie
        alphaMS(tag) = alpha_FMS + n_trop*tz;        %RA der mittl. Sonne
        MPG(tag) = 0;
        for nz=1:5                                  %MPG
            MPG(tag)  = MPG(tag)+ rad2deg(fac(nz)*sind(nz*M(tag)));
        end
        LS(tag)  = alphaMS(tag) + MPG(tag);         %Planetozentr. Laenge
        R2E(tag) = 0;
        for nz=1:10                                  %Reduktian auf Aequator
            faktor   = (-(tand(epsP/2))^2)^nz;
            R2E(tag) = R2E(tag)+rad2deg(faktor*sind(2*nz*LS(tag))/nz);
        end
        ZGL(PlNr, tag) = -4*(MPG(tag)+ R2E(tag));  %ZGL
        Dec(PlNr, tag) = asind(sind(epsP)*sind(LS(tag)))+ 0.25*sind(LS(tag));  %Dekl.
    end
end

%%
%%-----------------------------------------------------------------------%%
% Graphische Ausgabe
%

% Aufbereitung zum Ausgeben
for PlNr=PlSelect1:PlSelect2
    for m = NrDays(PlNr)+1:length(Dec(PlNr,:))
       ZGL(PlNr,m) = NaN;  
       Dec(PlNr,m) = NaN;
    end
end

% Bild Zeitgleichung
figure(); 
for PlNr =PlSelect1:PlSelect2 %Darstellung ab Erde bis Neptun
    subplot(2,(PlSelect2-2)/2,PlNr-2);
    plot(dtv(PlNr,:),ZGL(PlNr,:),...
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
    if PlNr == 8
        ylim([-40 40]);
    end
    if PlNr == 7
        ylim([-100 100]);
        xlim([-12 +12]);
        xlabel('ZGL in h');
    end
    set(gca,'Fontsize',16);
end



% Druckausgabe
fprintf('\n %s Keplerlösung \n', string(dt1) );
fprintf('Planet \n \n');
for PlNr =PlSelect1:PlSelect2
    %Länge eines Orbits/Jahres in Julianischen Tagen 
    n_trop= BaPadot.L0dot(PlNr);
    NrJDpOrbit = 365.25*36000/n_trop;           %Erdtage eines Jahres 

    %Länge eines Sonnentages (sol) und Länge eines Orbits/Jahres in sols
    SidTag        = BaPa.SidTagP(PlNr);         %Siderischer Tag
    SolTag        = SidTag/(1-SidTag/(NrJDpOrbit*86400)); 
                                                %Solartag auf Planeten(sol)
    NrSolDpOrbit  = NrJDpOrbit*86400/SolTag;    %Solartage eines Jahres 
    SolSidRatio   = SolTag/SidTag;              %Verhältnis Sol/Sid Tag Mars
    SolJDRatio    = SolTag/EarthDay;            %Verhältnis Sol/Erdtag
    SidJDRatio    = SidTag/EarthDay;            %Verhältnis Sid/Erdtag

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
% Functions
%
%--------------------------------------------------------------------------
% RotAxis: Berechnet die Orientierung des planetozentrischen Bezugssystems +
%         liefert die Rotationsrichtung
% Eingabe:
%   PlNr      Identifiziert den Planeten 
%   t         Zeit in Julianischen Jahrhunderten seit J2000
% Ausgabe:
%   alphaA,
%   deltaA    Richtunsgvekroren der Polachse des Planeten
%   Sense     Rotationsrichtung (direkt oder retrograd)
%--------------------------------------------------------------------------
function [alphaA, deltaA, Sense] = RotAxis (PlNr,t)
  % Rektaszension und Deklination der Rotationsachse bezogen auf
  % Aequator und Ekliptik zur Epoche J2000; Orientierung des
  % Nullmeridians
  % Position der Planetenachse in geozentrisch äquatorialen Koord. Ref. [3]

  switch (PlNr)
    case 1
        RA  = 281.01  -   0.033*t;  Dec =  61.45  -   0.005*t;    
    case 2
        RA  = 272.76; Dec =  67.16; 
    case 3
        RA  =   0.00  -   0.641*t;  Dec =  90.00  -   0.557*t;
    case 4
        RA  = 317.681 -   0.108*t; Dec =  52.886 -   0.061*t;
    case 5
        RA  = 268.05  -   0.009*t; Dec =  64.49  +   0.003*t;
    case 6
        RA  =  40.589 -   0.036*t; Dec =  83.537 -   0.004*t;  
    case 7
        RA  = 257.311; Dec = -15.175; 
    case 8
        N   = deg2rad(357.85+52.316*t);    
        RA  = 299.36  + 0.70*sin(N); 
        Dec =  43.46  - 0.51*cos(N); 
    case 9
        RA  = 313.02; Dec =   9.09;
  end
  alphaA  = RA; 
  deltaA  = Dec; 
  % Rotationsrichtung
  if ( PlNr==2 || PlNr==7 || PlNr==10 )
    Sense = +1;
  else
    Sense = -1;
  end
end 
        

