% -------------------------------------------------------------------------
% Sonnenfinsternis03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
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
% Berechnet die lokalen Umstaende einer Sonnenfinsternis.
% Fuer die Koordinaten von Sonne und Mond werden genaue Reihen verwendet.
% Input ist die ungefaehre Finsterniszeit. 
%
% Es wird eine geographische Position erwartet, die im Finsternisbereich
% liegt. Evtentuell dazu SolarEclips02 verwenden.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
SunData = PerturbImport('SonnePos.csv');

% Startdaten
% z.B.
% Oberkochen 
phi  =  48.1;
lambda = 10.9;
d_0   = "11-Aug-1999 09:30:00"; %(UT)

% oder
% Snake River 2017
% phi  =  44.3 ;
% lambda = -117.1;
% d_0   = "21-Aug-2017 17:30:00"; %(UT)

SITE = Site(phi);
dt_0  = datetime(d_0,'InputFormat','dd-MMM-yyyy HH:mm:ss');
T_0   = juliandate(dt_0);
DelTSec = ETminusUT(T_0);
t_F = Jd2JJht(T_0); %Zeit in Julianischen Jahrhunderten seit JD2000
shadow = 'partiell';

% Eigentlich Suchprogramm nach den Kontaktzeiten
CONTACTS = contacts(t_F, DelTSec, lambda, SITE, SunData);

if ~strcmpi(CONTACTS.PHASE,'keine Finsternis')
    outstr(:) = string(4);
    for i = 1:4 
           if isnan(CONTACTS.times(i))
               outstr(i,:) = '----- ';
           else
               zeit(i) = datetime(JJht2Jd(CONTACTS.times(i)), ...
                               'ConvertFrom','juliandate');
               outstr(i,:) = string(zeit(i),'HH:mm:ss')+" UT";
           end
end
    zeit =   datetime(JJht2Jd(CONTACTS.t_MAX),'ConvertFrom','juliandate');
    outstr_max = string(zeit,'dd-MMM-yyyy HH:mm:ss')+" UT";
    fprintf('\n Lokale Umstaende der Sonnenfinsternis fuer  \n ');
    fprintf('\n Datum: %s   Breite: %s° Laenge: %s°  \n ',...
              string(dt_0, 'dd-MMM-yyyy'), ...
              num2str(phi, '%+5.2f'), num2str(lambda,'%+5.2f'));
    fprintf('Delta T = ET - UT : %6.3f \n \n ', DelTSec );
    fprintf('Finsternis : %s \n ', CONTACTS.PHASE);
    fprintf('\n Maximum der Verfinsterung: %s  \n ',outstr_max);
    for i = 1:4 
        fprintf('%1d. Kontakt %s   \n ',i, outstr(i,:));
    end
    fprintf('\n Groesse Finsternis (bedeckter Durchmesserteil): %s \n',...
             num2str(CONTACTS.MAG));
    fprintf(' Verfinsterungsgrad [1 = Flaeche der Sonne]    : %s \n',...
             num2str(CONTACTS.OBSCUR));
end
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%__________________________________________________________________________
%%
%  Funktionen

%--------------------------------------------------------------------------                                                                           
% CONTACTS                                                                                                                                               
function CONTACTS = contacts(t_F, DelTSec, lambda, SITE, SunData)
%
%   Bestimmt ausgehend vom Neumondzeitpunkt fuer einen bestimmten             
%   Beobachtungsort die Phase, Groesse und Kontaktzeiten einer                
%   Sonnenfinsternis.                                                         
%
% Input:
%   t_F                Ungefaehre Neumondzeit in UT                           
%   DelTSec            Differenz Ephemeridenzeit-Weltzeit (Sekunden)    
%   lambda,SITE        Geozentrische Koordinaten des Beobachters        
%   SunData            Parameter Sonnenkoordinaten                                  
% Output:
%   CONTACTS.times     Zeiten des 1. bis 4. Kontakts (UT)               
%   CONTACTS.tmax      Zeit der maximalen Verfinsterung (UT)            
%   CONTACTS.MAG       Groesse der Finsternis (bedeckter Bruchteil      
%                      des Sonnendurchmessers im Maximum)               
%   CONTACTS.OBSCUR    Verfinsterungsgrad [1 = Flaeche der Sonne]       
%   CONTACTS.PHASE     Phase im Maximum der Finsternis                  
%                                                                             
% Alle Zeiten in julianischen Jahrhunderten seit J2000.                       
% Beschreibung:                                                               
%                                                                             
% Zur Bestimmung der Kontaktzeiten werden die Nullstellen der Funktion        
% f(t)=D(t)**2-L(t)**2 bestimmt, wobei D den Abstand des Beobachters von      
% der Schattenachse und L den Radius des Schattenkegels bezeichnet.           
% Fuer f(t)=0 liegt der Beobachter genau auf dem Rand des Schattenkegels      
% und sieht dabei, wie sich die Raender von Sonne und Mond scheinbar          
% beruehren.                                                                  
%                                                                             
% Mit Hilfe der quadratischen Interpolation wird zunaechst nach dem 1. und    
% 4. Kontakt gesucht, bei dem der Beobachter auf dem Halbschattenkegel        
% liegt. Gleichzeitig kann damit der Zeitpunkt der maximalen Finsternis       
% festgestellt werden, bei der der Abstand zwischen Beobachter und            
% Schattenachse minimal wird.                                                 
%                                                                             
% Ist die Finsternis wenigstens partiell, so berechnet CONTACTS fuer den      
% Zeitpunkt des Finsternismaximums die Groesse der Finsternis und die         
% Bedeckung der Sonne durch den Mond.                                         
%                                                                             
% Ergibt sich aus dem Abstand des Beobachters von der Schattenachse und       
% dem Kernschattendurchmesser zum Zeitpunkt des Maximums, dass die            
% Finsternis total oder ringfoermig ist, so koennen anschliessend die         
% Zeiten des 2. und 3. Kontakts bestimmt werden. Da die Dauer der totalen     
% oder ringfoermigen Phase nie laenger als ca.13 Minuten dauert, liegen der   
% 2. und 3. Kontakt jeweils in einem entsprechenden Intervall vor bzw. nach   
% dem Maximum.    
%                                                                             
%-------------------------------------------------------------------------- 

    shadow = 'partiell';
    t_UT = t_F;
    myfunc = @(x, a, b, c, d, e) sha_dist(shadow, t_UT, DelTSec,...
                                             lambda, SITE, SunData);
    a = shadow;         % parameter1
    b = DelTSec;        % parameter2
    c = lambda;         % parameter3
    d = SITE;           % parameter4
    e = SunData;        % parameter5
    func = @(x)myfunc(a, x, b, c, d, e);  % function of x alone

    RANGE = 5.0;          % Suchzeitraum = +/-(RANGE+DTAB) [h]            
    DTAB  = 1/120;        % Schrittweite [h]                              
    CENT  = 876600.0;     % 24*36525 Stunden pro Jahrhundert              
    EPS   = 1.0E-10;      % Toleranz fuer Kontaktzeitensuche (ca. 0.3s) 
    t_Max = NaN;
    % Initialisierung   
    for i = 1:4 
        CONTACTS.times(i) = t_F;
    end
    %----------------------------------------------------------------------  
    % Suche nach dem 1. und 4. Kontakt mit quadratischer Interpolation.       
    %                                                                         
    % Der Verlauf des Halbschattenabstandes zum Beobachter in der Funda-      
    % mentalebene wird durch je drei aequidistante Funktionswerte             
    % PD_MINUS, PD_0 und PD_PLUS im zeitlichen Abstand DTAB dargestellt.      
    % Durch diese drei Punkte wird eine interpolierende Parabel gelegt,       
    % deren Nullstellen bestimmt werden. Dieses Verfahren liefert mit         
    % geringem Aufwand die Zeiten des 1. und 4. Kontakts, weil die inter-     
    % polierte Funktion (das Quadrat des Halbschattenabstand) ausreichend     
    % glatt ist.                                                              
    %----------------------------------------------------------------------  
    N_CONTACTS = 0;
    t_UT    = t_F + (-RANGE-DTAB)/CENT;
    PD_PLUS = sha_dist('partiell', t_UT, DelTSec,lambda, SITE, SunData);
    t_UT    = t_F + (-RANGE-2.0*DTAB)/CENT;
    t_MAX   = 0;
    index = 0;
    while (t_UT < t_F+ RANGE/CENT) && (N_CONTACTS < 2)
        index =index+1;
        % Naechsten Zeitpunkt berechnen   
        t_UT = t_UT + 2.0*DTAB/CENT;
        % Quadrat des Schattenabstandes zu den Zeiten T_UT-DTAB, T_UT   
        % und T_UT+DTAB berechnen und interpolieren                     
        PD_MINUS = PD_PLUS;
        PD_0     = sha_dist ('partiell', t_UT, DelTSec, lambda, ...
                   SITE, SunData);
        PD_PLUS  = sha_dist ('partiell', t_UT + DTAB/CENT, DelTSec,...
                   lambda, SITE, SunData);
        %Quadratische Interpolation fuer Nullstellen(normiert auf [-1,1])
        Roots = quad(PD_MINUS,PD_0,PD_PLUS);
        
        %Eventuelle Zwischenausgabe
        
        T_UT = JJht2Jd(t_UT);
        zeit = datetime(T_UT,'ConvertFrom','juliandate');  
        temp1(index)=zeit;
        outstr = string(zeit,'dd-MMM-yyyy HH:mm')+" UT";
        
        T_MAX =JJht2Jd(t_MAX);
        zeit =   datetime(T_MAX,'ConvertFrom','juliandate');
        outstr1 = string(zeit,'dd-MMM-yyyy HH:mm')+" UT";
        % Anzahl der bereits gefundenen Kontakte speichern und   
        % Kontaktzeiten berechnen    
        N_ROOTS = Roots.NZ;
        XE = Roots.XE;
        YE = Roots.YE;
        ROOT1 = Roots.ZERO1;
        ROOT2 = Roots.ZERO2;
        N_CONTACTS = N_CONTACTS+N_ROOTS;
        switch (N_ROOTS)
            case 1 
                if  (N_CONTACTS==1) %1. Kontakt   
                    CONTACTS.times(1) = t_UT+(DTAB*ROOT1)/CENT;
                else
                    CONTACTS.times(4) = t_UT+(DTAB*ROOT1)/CENT;% 4. Kontakt   
                end
            case 2
                CONTACTS.times(1) = t_UT+(DTAB*ROOT1)/CENT;  % 1. Kontakt   
                CONTACTS.times(4) = t_UT+(DTAB*ROOT2)/CENT;  % 4. Kontakt   
        end
        % Falls das Minimum des interpolierten Schattenabstandes im Such-       
        % intervall liegt, wird aus seiner Lage der Zeitpunkt der maximalen     
        % Verfinsterung ausreichend genau, ohne weitere Iteration bestimmt.  
        if (-1.0<XE) && (XE<1.0) 
            t_MAX = t_UT + (DTAB*XE)/CENT;
        end
    end     
    % Finsternistyp festlegen
    if N_CONTACTS == 0
            PHASE = 'keine Finsternis'; % Keine Finsternis gefunden            
    else
            PHASE = 'partiell'; % Finsternis ist mindestens partiell
    end 
         
    %----------------------------------------------------------------------
    %                                                                         
    % Verfinsterungsgrad und Groesse der Finsternis                           
    %                                                                         
    %----------------------------------------------------------------------
    OBSCUR = 0;
    MAG =0;
    if strcmp(PHASE,'partiell')
        % Koordinaten von Schatten und Beobachter auf der Fundamentalebene   
        % im Maximum der Finsternis   
        BESSEL = bessel(t_MAX, DelTSec, SunData);
        Obs = Observer(t_MAX, lambda, SITE, BESSEL);
        % Abstand des Beobachters von der Schattenachse   
        M = sqrt((Obs.xi-BESSEL.SH(1))^2 + (Obs.eta-BESSEL.SH(2))^2);
        % Penumbra- und Umbraradius am Ort des Beobachters    
        % (LL2<0 fuer totale Finsternis!)   
        LL1 = BESSEL.l1 - Obs.zeta*tan(BESSEL.f1);  
        LL2 = BESSEL.l2 - Obs.zeta*tan(BESSEL.f2);  
        % Art der Finsternis   
        if (M < LL2)
            PHASE = 'ringfoermig';
        end
        if M < - LL2 
               PHASE = 'total';
        end
        % Groesse der Finsternis   
        if strcmp(PHASE,'partiell')
            MAG = (LL1-M)/(LL1+LL2);    % Im HS  
        else
            MAG = (LL1-LL2)/(LL1+LL2);  % Im KS      
        end
    end
    %Verfinsterungsgrad berechnen   
    switch PHASE
         case 'keine Finsternis' 
            OBSCUR = 0.0;
         case 'partiell' 
            B =  acos((LL1*LL2+M*M)/(M*(LL1+LL2)));
            C =  acos((LL1*LL1+LL2*LL2-2*M*M)/(LL1*LL1-LL2*LL2));
            A =  pi-(B+C);
            S = (LL1-LL2)/(LL1+LL2);
            OBSCUR = (S*S*A + B - S*sin(C))/pi;
         case 'ringfoermig'
            S = (LL1-LL2)/(LL1+LL2);
            OBSCUR = S*S;
         case 'total'
            OBSCUR = 1.0;
    end
    
    %----------------------------------------------------------------------  
    % Die Zeiten des 2. und 3. Kontaktes ausgehend vom Zeitpunkt der maxi-    
    % malen Verfinsterung mittels Pegasus-Verfahren suchen. Suchintervall     
    % [T_MAX-DT,T_MAX] fuer 2. Kontakt und [T_MAX,T_MAX+DT] fuer 3.Kontakt.   
    % Die Fehlertoleranzgrenze EPS entspricht ca. 0.3 Sekunden.               
    %---------------------------------------------------------------------- 
    if strcmp(PHASE,'ringfoermig') || strcmp(PHASE,"total")
        
       DT0 = 0.10/CENT; % +- 6 min in julianischen Jahrhunderten   
       DT = 0.10/60/12/CENT; % +- 5 sec in julianischen Jahrhunderten   
       shadow = 'total';
       ts = t_MAX-DT0;
       y1 = sha_dist(shadow,ts , DelTSec, lambda, SITE, SunData);
       y2 = sha_dist(shadow,ts+DT, DelTSec, lambda, SITE, SunData);
       ncontacts34 =0;
       zaehler = 1;
       while ncontacts34 <2 && zaehler < 1000
           if y1*y2 < 0 
               if y1 > 0
                   CONTACTS.times(2) = ts + DT/2;
                   ncontacts34 = ncontacts34 +1;
               else
                   CONTACTS.times(3) = ts + DT/2;
                   ncontacts34 = ncontacts34 +1;
               end
           end  
           y1 = y2;
           ts = ts + DT;
           y2 = sha_dist(shadow,ts+DT, DelTSec, lambda, SITE, SunData);
           zaehler = zaehler +1;
%            fprintf( ' %s  : %f \n ', string(datetime(JJht2Jd(ts),...
%                     'ConvertFrom','JulianDate'),'HH:mm:ss'), y1);
       end
    else
       CONTACTS.times(2) = NaN;
       CONTACTS.times(3) = NaN;
    end
 
    % Rueckgabewerte
    CONTACTS.OBSCUR = OBSCUR;
    CONTACTS.t_MAX = t_MAX;
    CONTACTS.MAG = MAG;
    CONTACTS.PHASE = PHASE;
end
%--------------------------------------------------------------------------
% Ende Funktion contacts
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
function SITE = Site(PHI)
%
%- Output: SITE:  berechnet geozentr. aus geograph. Beobachterkoordinaten     
%-          Site.c:  r * cos(phi') (geozentrisch; in Erdradien)                   
%-          Site.s:  r * sin(phi') (geozentrisch; in Erdradien)                   
%- Input:  PHI:   geographische Breite (in Grad)                               
%--------------------------------------------------------------------------  
    ex2 = 0.006694;        % ex^2=f(2-f) mit Erdabplattung f=1/298.257 
    sinphi = sind(PHI);   
    cosphi = cosd(PHI);
    N = 1.0/sqrt(1.0-ex2*sinphi*sinphi);
    SITE.c = N*cosphi;  
    SITE.s = (1.0-ex2)*N*sinphi;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Ende Funktion Site
%--------------------------------------------------------------------------


%%
function Roots = quad(Y_MINUS,Y_0,Y_PLUS)
%--------------------------------------------------------------------------  
% QUAD: 
%     bestimmt die Parabel durch 3 Punkte _MINUS), (0,Y_0) und (1,Y_PLUS),                                 
%     die nicht auf einer Geraden liegen.                                
%     Y_MINUS,Y_0,Y_PLUS: drei y-Werte                                       
%     XE,YE   : x und y Wert des Extremums der Parabel                       
%     ZERO1   : erste Nullstelle im Intervall [-1,+1] (fuer NZ=1,2)          
%     ZERO2   : zweite Nullstelle im Intervall [-1,+1] (fuer NZ=2)           
%     NZ      : Zahl der Nullstellen der Parabel im Intervall [-1,+1]        
%--------------------------------------------------------------------------
    NZ = 0;
    A  = 0.5*(Y_MINUS+Y_PLUS)-Y_0; 
    B  = 0.5*(Y_PLUS-Y_MINUS); 
    C  = Y_0;
    ZERO1 = 0;
    ZERO2 = 0;
    XE = -B/(2.0*A); 
    YE = (A*XE + B) * XE + C;
    DIS = B*B - 4.0*A*C;  % Diskriminante von y = axx+bx+c 
    if (DIS >= 0)         % Parabel hat Nullstellen        
        DX    = 0.5*sqrt(DIS)/abs(A); 
        ZERO1 = XE-DX; 
        ZERO2 = XE+DX;
        if abs(ZERO1) <= +1.0 
            NZ = NZ + 1;  
        end
        if abs(ZERO2) <= +1.0
            NZ = NZ + 1;
        end
        if  ZERO1<-1.0
            ZERO1 = ZERO2;
        end
    end
    Roots.XE = XE;
    Roots.YE = YE;
    Roots.ZERO1 = ZERO1;
    Roots.ZERO2 = ZERO2;
    Roots.NZ    = NZ;
end
%--------------------------------------------------------------------------
% Ende Funktion quad
%--------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
%
% SHA_DIST = sha_dist(shadow, t_UT,DelTSec, lambda, SITE, SunData)
%
%  Berechnet die Schattenabstandsfunktion f(t)=D(t)**2-L(t)**2, wobei D  
%  den Abstand des Beobachters von der Schattenachse und L den Radius des
%  Schattenkegels bezeichnet. Fuer f(t)=0 liegt der Beobachter genau auf
%  dem Rand des Schattenkegels. (Einheit Erdradien**2).
%
% Eingabe:
%  t_UT        Berechnungszeit in Julianischen Jahrhunderten UT seit J2000
%  shadow      Umbra oder Penumbra als string "total" oder "partiell"
%  DelTSec     Differenz zwischen Ephemeridenzeit und Weltzeit in [s]
%  SITE        Geozentrische Beobachterkoordinaten
%  SunData     Parameter fuer Sonnenkoordinaten-Berechnung
% Ausgabe:
%  SHA_DIST    Abstand Beobachter-Schattenmitte in [Erdradien]
%
function SHA_DIST = sha_dist(shadow, t_UT, DelTSec, lambda, SITE, SunData)
    K_Moon = 0.2725076; %RM/RE
    % Beobachter- und Schattenkoordinaten im System der Fundamentalebene   
    BESSEL = bessel  (t_UT, DelTSec, SunData);
    Obs    = Observer(t_UT, lambda , SITE, BESSEL);
    % Schattenradius am Ort des Beobachters (LL<0 fuer totale Finsternis)   
    if strcmp(shadow,'partiell')
         LL = BESSEL.l1-Obs.zeta*tan(BESSEL.f1);  
         SHA_DIST = 1*((BESSEL.SH(1)-Obs.xi)^2 + ...
                       (BESSEL.SH(2)-Obs.eta)^2) - 1* LL*LL;
    else
         LL = BESSEL.l2-Obs.zeta*tan(BESSEL.f2);  
         SHA_DIST = 1*((BESSEL.SH(1)-Obs.xi)^2 + ...
                       (BESSEL.SH(2)-Obs.eta)^2) - 1* LL*LL;

    end 
end
%-------------------------------------------------------------------------
% Ende Funktion sha_dist
%-------------------------------------------------------------------------

%%
%--------------------------------------------------------------------------
% Bessel: Berechnet die Orientierung der Fundamentalebene, die Koordinaten
%         der Schattenmitte und die Parameter des Schattenkegels
%
% Eingabe:
%   t_UT        Berechnungszeit in Julianischen Jahrhunderten UT seit J2000
%   DelTSec     Differenz zwischen Ephemeridenzeit und Weltzeit in [s]
%   SDL, SDR    Sonnenkoordinaten 
% Ausgabe:
%   BESSEL      
%   Struktur mit:
%   IJK         Basisvektoren der Fundamentalebene (aequatorial)
%   XSh, YSh    Koordinaten der Schattenmitte auf der Fundamentalebene 
%               in [Erdradien]
%   F1          Halber Oeffnungswinkel des Penumbra-Schattenkegels in [rad]
%   L1          Radius des Penumbra-Schattenkegels auf der Fundamentalebene
%               in [Erdradien]
%   F2          Halber Oeffnungswinkel des Umbra-Schattenkegels in [rad]
%   L2          Radius des Umbra-Schattenkegels auf der Fundamentalebene
%               in [Erdradien]
%--------------------------------------------------------------------------
function BESSEL = bessel(t_UT, DelTSec, SunData)

% Berechnet Orientierung (IJK) der Hauptebene, Koordinaten Schattenmitte
% Vec SH und Parameter des Schattenkegels W1, D1, W2, D2
    K_Moon = 0.2725076; %RM/RE
    K_Sun  = 109.1227;  %RS/RE  
    tau = 8.32/1440/36525; %Lichtlaufzeit;
    epsE=deg2rad(23.43929111);


    t_ET = t_UT + DelTSec/86400/36525;
    SunPos = SonneExakt(t_ET-tau,SunData,epsE,'RE'); %Beachte tau ! 
    MoonPos = MondExakt(t_ET,epsE,'RE');
    
    X_M = MoonPos.geo; %kartes. aequatoriale Koord. v. Mond und Sonne
    X_S = SunPos.geo;
    
    X_MS   = X_S-X_M;
    R_MS   = norm(X_MS);
    E_XS   = X_MS/R_MS; %Richtungsvektor Sonne-Mond
    %Neues Koordinatensystem IJK fuer Fundamentalebene
    BESSEL.IJK3   = E_XS;
    DIST   = sqrt(E_XS(1)*E_XS(1)+E_XS(2)*E_XS(2));
    BESSEL.IJK1   = [-E_XS(2)/DIST;+E_XS(1)/DIST;0];
    BESSEL.IJK2   = cross(BESSEL.IJK3, BESSEL.IJK1);
    %Mondkoordinaten auf der Fundamentalebene
    BESSEL.SH    = [dot(X_M, BESSEL.IJK1);dot(X_M, BESSEL.IJK2); ...
                    dot(X_M, BESSEL.IJK3)];
    %Entfernung Sonne-Mond
    %Schattenkegel
    BESSEL.f1    = asin((K_Sun+K_Moon)/R_MS);
    BESSEL.f2    = asin((K_Sun-K_Moon)/R_MS);

    BESSEL.l1    = BESSEL.SH(3)*(K_Sun+K_Moon)/R_MS + K_Moon;
    BESSEL.l2    = BESSEL.SH(3)*(K_Sun-K_Moon)/R_MS - K_Moon;
    
    BESSEL.v1    = K_Moon*DIST/(K_Sun+K_Moon);
    BESSEL.v2    = K_Moon*DIST/(K_Sun-K_Moon);
end
%-------------------------------------------------------------------------
% Ende Funktion bessel
%-------------------------------------------------------------------------

%--------------------------------------------------------------------------
function OBS = Observer(t_UT, lambda, SITE, BESSEL)
%
% Observer: 
% Projektion des Beobachtungsortes auf die Fundamentalebene
%
% Eingabe:
%   t_UT        Berechnungszeit in Julianischen Jahrhunderten UT seit J2000
%   Bessel.IJK  Basisvektoren der Fundamentalebene (aequatorial)
%   SITE        Geozentrische Beobachterkoordinaten
% Ausgabe:      Koordinaten des Beobachters auf der Fundamentalebene 
%               in [Erdradien]
%--------------------------------------------------------------------------
    % Stundenwinkel des Beobachters in Grad
    Tau = 15*LMST(t_UT*36525+51544.5,lambda);
    % Äquatoriale kartesische Korodinaten des Beobachters
    OBS.equ = [SITE.c*cosd(Tau);SITE.c*sind(Tau);SITE.s];
    % Projektion in die Hauptebene
    OBS.xi   = dot(OBS.equ, BESSEL.IJK1);
    OBS.eta  = dot(OBS.equ, BESSEL.IJK2);
    OBS.zeta = dot(OBS.equ, BESSEL.IJK3);  
end
%--------------------------------------------------------------------------


%%
%_________________________________________________________________________
%
% Ende Programm
%_________________________________________________________________________




