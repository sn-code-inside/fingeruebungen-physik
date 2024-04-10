% -------------------------------------------------------------------------
% SonnenaufgangBaabe.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Sonnenaufgang Baabe 
% Berechnung von Alt und Az bei SA und SU 
% Programm berechnet die Alt_Az_Kurve beim Sonnenaufgang 
% auf Basis einer (vereinfachten) MPG- Keplerlösung 
% für Baabe Frühlingsanfang und zur Sommersonnenwende.
% Dabei wird einmal die "astrometrisch=wahre" Höhe und einmal die durch
% Refraktion u.a. Effekte bedingte "scheinbare" Höhe berechnet.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

pi2 = 2*pi;                   % 2*pi

% Koordinaten Baabe/Ostsee
lambda=13.7;
phi=54.36;
Zone=0;
phistr    = sprintf('%+4.2f ',phi);
lambdastr = sprintf('%+4.1f',lambda);

phi=deg2rad(phi);
sinphi=sin(phi);
cosphi=cos(phi);

hr=-(50/60);
ds=-16/60;
hrn=hr+8/3600;

dt1 = datetime('2000-01-01 12:00:00');
T0 = juliandate(dt1) - Zone/24; % Julianisches Datum
MJuDa0 = juliandate(dt1,'modifiedjuliandate')- Zone/24;% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
MJuDa0 = MJuDa0-0.5;

Alt1 = zeros(5,120);
Alt1k= zeros(5,120);
Alt2 = zeros(8,5,120);

sunrise1 = zeros(5);
sunset1  = zeros(5);
nullinie = zeros(120);
timeshift = zeros(8,5);


tag(2)=datetime(2000,3,20);
tag(3)=datetime(2000,6,21);

titlestr   =  strings([5,25]);
titlestr(2,:)= string(tag(2),'dd.MM.');
titlestr(3,:)= string(tag(3),'dd.MM.');

%_________________________________________________________________________
% Hauptprogramm

% Schleife ueber 2 Zeitpunkte
  
for n=2:3
  MJuDa = MJuDa0 + day(tag(n),'dayofyear');
  T=T0+day(tag(n),'dayofyear');
  hr=deg2rad(-(50/60));
% Berechnung SA nach MPG-Keplerlösung mit quadratischer Interpolation
  [sunrise1(n), sunset1(n), xRise(n), xSet(n), xAbove(n)] = FindRiseSet(MJuDa,lambda,phi);
  sunrise1(n) = (wrapTo360((sunrise1(n)+Zone)*15)/15);
  if xRise(n)
        MJuDaRise=MJuDa+sunrise1(n)/24;
        zeit0 = -15/60-1/3600;
        zeit = zeros(120);
        for m= 1:120
           zeit(m)    = zeit0 + 15*m/3600;
           Alt1(n,m)  = asind(SinAlt(MJuDaRise,zeit(m),lambda,cosphi,sinphi));
           Az1(n,m)   = wrapTo180(atan2d(TanAz(MJuDaRise,zeit(m),lambda,cosphi,sinphi),1));
           if Az1(n,m)< 0 
               Az1(n,m) = Az1(n,m)+180;
           end
           has=Alt1(n,m);
           if has > hrn
              refraction=(ds+1.02/tand(has+10.3/(has+5.11))/60); 
           else
              refraction=(ds+1.02/tand(hrn+10.3/(hrn+5.11))/60); 
           end
           refraction=-hrn;
           Alt1k(n,m) = Alt1(n,m)+refraction;
        end
  else
       for m= 1:120
           Az1(n,m)   = pi;
           Alt1(n,m)  = 0.5*xAbove(n);
           Alt1k(n,m) = Alt1(n,m);
       end
  end
  Azmin(n)=round(min(Az1(n,:)));
end

%------------------------------------------------------------------------------
%  Ausgabe
 

header1='Winkel bei Sonnenaufgang'; 
figure('Name',header1);
yaxsc=60;
for k=2:3
    subplot(1,2,k-1);
    plot(Az1(k,:), Alt1(k,:)*yaxsc,Az1(k,:), Alt1k(k,:)*yaxsc,Az1(k,:), nullinie,'g:','LineWidth',2);
    xlim([Azmin(k)  Azmin(k)+6]);
    ylim([-1*yaxsc 1*yaxsc]);
    title(titlestr(k)+' {\phi} = '+phistr+'°'+' {\lambda} = '+lambdastr+'°');
    grid on;
    grid minor;
    xlabel('Azimuth in °');
    ylabel('Höhe in arcmin');
    legend('astrometrisch','scheinbar','location','northwest');
    legend boxoff;
end

header1='Zeitdelay zw. wahrem und scheinbarem SA (- = früher)'; 
figure('Name',header1);
for k=2:3
    subplot(1,2,k-1);
    plot(Alt1(k,:),zeit*60,Alt1k(k,:),zeit*60,'LineWidth',2);
    hold on;
    plot(nullinie, zeit*60,'Color','green','LineStyle',':','LineWidth',2);
    xlim([-1 1]);
    ylim([-2 10]);
    title(titlestr(k)+' {\phi} = '+phistr+'°'+' {\lambda} = '+lambdastr+'°');
    grid on;
    grid minor;
    ylabel('Minuten');
    xlabel('Höhe in °');
    legend('astrometrisch','scheinbar','location','southeast');
    legend boxoff;
end


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%Hilfsfunktion zur Berechnung des Azimuth
function tanaz=TanAz(MJT, hour, lambda, cosphi, sinphi)
  MJD=MJT+hour/24;
  T=Jd(MJD);
  eps0=deg2rad(23.43929111);      % Schiefe der Ekliptik
  [xRA,xDec]=KeplerSonne(T,eps0);
  tau=GMST(MJD)+deg2rad(lambda)-xRA;
  tau   = wrapToPi(tau);
  tanaz = cos(xDec)*sin(tau)/(sinphi*cos(xDec)*cos(tau)-cosphi*sin(xDec))+0.000001;
end

