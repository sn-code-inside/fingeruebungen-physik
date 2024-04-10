% -------------------------------------------------------------------------
% SonnenaufgangZeit.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Sonnenaufgangsberechnung Alt Az 
% Programm berechnet die Alt_Az_Kurve beim Sonnenaufgang
% auf Basis einer (vereinfachten) MPG-Keplerlösung 
% als Funktion verschiedener geographischer Breiten und für
% das Äquinoktium im Frühjahr und die Sommersonnwende.
% Dabei wird  die "astrometrisch" wahre Höhe und die durch
% Refraktion u.a. Effekte "scheinbare" beobachtete Höhe berechnet.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
pi2 = 2*pi;                   

lambda=15;
lambdastr = sprintf('%+4.1f',lambda);
Zone=0;

hr=-(50/60);
ds=-16/60;
hrn=hr+8/3600;

dt1 = datetime('2000-01-01 12:00:00');
T0 = juliandate(dt1) - Zone/24; % Julianisches Datum
MJuDa0 = juliandate(dt1,'modifiedjuliandate')- Zone/24;% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
MJuDa0 = MJuDa0-0.5;


%------------------------------------------------------------------------------

Alt1 = zeros(6,120);
Alt1k= zeros(6,120);
Az1  = zeros(6,120);  
sunrise1 = zeros(6);
sunset1 = zeros(6);
nullinie = zeros(120);

tag(1)=datetime(2000,3,20);  %Frühlingsequinox
tag(2)=datetime(2000,6,21);  %Sommersonnenwende

titlestr   =  strings([2,25]);
titlestr(1,:)= string(tag(1),'dd.MM.');
titlestr(2,:)= string(tag(2),'dd.MM.');

phistr   =  strings; %verschiedene Breitengrade
for L=1:3
    phi(L) = 5+ (L-1)*30;
    phistr(L,:) = sprintf(' %+4.2f',phi(L));
    phirad(L) = deg2rad(phi(L));
end

%_________________________________________________________________________
% Hauptprogramm


% Schleife ueber 2 Zeitpunkte und 3 geographische Breiten
for L=1:3
    sinphi=sin(phirad(L));
    cosphi=cos(phirad(L));
    for n=1:2
        MJuDa = MJuDa0 + day(tag(n),'dayofyear');
        T=T0+day(tag(n),'dayofyear');
        % Berechnung SA nach MPG-Kepler mit quadratischer Interpolation
        [sunrise1(n), sunset1(n), xRise(n), xSet(n), xAbove(n)] = FindRiseSet(MJuDa,lambda,phirad(L));
        sunrise1(n) = (wrapTo360((sunrise1(n)+Zone)*15)/15);
        if xRise(n)
            MJuDaRise=MJuDa+sunrise1(n)/24;
            zeit0 = -15/60-1/3600;
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
                Alt1k(n,m) = Alt1(n,m)+refraction; %Berücksichtigung Refraktion
            end
        else
            for m= 1:120
                Az1(n,m)   = pi;
                Alt1(n,m)  = 0.5*xAbove(n);
                Alt1k(n,m) = Alt1(n,m);
            end
        end    
        Azmin(L,n)=round(min(Az1(n,:)));
    end
    outy1(L).tracey1=Alt1(:,:);
    outy2(L).tracey2=Alt1k(:,:);
    outx(L).tracex=Az1(:,:);
end

%------------------------------------------------------------------------------
%  Ausgabe
 
header1='Winkel bei Sonnenaufgang'; 
figure('Name',header1);
yaxsc=60;
for L=1:3
     for k=1:2
        subplot(2,3,(k-1)*3+L);
        plot(outx(L).tracex(k,:), outy1(L).tracey1(k,:)*yaxsc, outx(L).tracex(k,:),outy2(L).tracey2(k,:)*yaxsc,'LineWidth',2);
        hold on;
        x=[Azmin(L,k) Azmin(L,k)+6];
        y=[0 0];
        line(x,y,'Color','green','LineStyle','--','LineWidth',2)
        xlim([Azmin(L,k)  Azmin(L,k)+6]);
        ylim([-1*yaxsc 1*yaxsc]);
        title(titlestr(k)+' {\phi} = '+phistr(L,:)+'°'+' {\lambda} = '+lambdastr+'°');
        grid on;
        grid minor;
        xlabel('Azimut in °');
        if yaxsc==60 
            ylabel('Höhe in arc min');
        else
            ylabel('Höhe in °');
        end
        legend('astrometrisch','scheinbar','location','southeast');
        legend boxoff;
     end
end

header1='Zeitdelay zw. wahrem und scheinbarem SA (- = früher)';
figure('Name',header1);
for L=1:3
    for k=1:2
        subplot(2,3,(k-1)*3+L);
        plot(outy1(L).tracey1(k,:),zeit*60,outy2(L).tracey2(k,:),zeit*60,'LineWidth',2);
        hold on;
        plot(nullinie, zeit*60,'Color','green','LineStyle','--','LineWidth',2);
        xlim([-1 1]);
        ylim([-5 15]);
        title(titlestr(k)+' {\phi} = '+phistr(L,:)+'°'+' {\lambda} = '+lambdastr+'°');
        grid on;
        grid minor;
        ylabel('Minuten');
        xlabel('Höhe in °');
        legend('astrometrisch','scheinbar','location','southeast');
        legend boxoff;
    end
end


% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%Hilfsfunktion zur Berechnung des Azimut
function tanaz=TanAz(MJT, hour, lambda, cosphi, sinphi)
  MJD=MJT+hour/24;
  T=MJD+2400000.5;
  eps0=deg2rad(23.43929111);      % Schiefe der Ekliptik
  [xRA,xDec]=KeplerSonne(T,eps0);
  tau=GMST(MJD)+deg2rad(lambda)-xRA;
  tau   = wrapToPi(tau);
  tanaz = cos(xDec)*sin(tau)/(sinphi*cos(xDec)*cos(tau)-cosphi*sin(xDec))+0.000001;
end


function sinusalt=SinAlt(MJT, hour, lambda, cosphi, sinphi)
% Eingabe: 
%   MJT     = modifiziertes julianisches Datum 
%   hour    = Zeitpunkt in Stunden 
%   lambda  = geogr. Länge 
%   cosphi  = Cosinus der geogr. Breite
%   sinphi  = Sinus der geogr. Breite
% Ausgabe: 
%   sinusalt = Sinus der Sonnenhöhe 

  MJD=0;
  MJD=MJT+hour/24;              % Umrechnung der Zeit 
  eps=deg2rad(23.43929111);     % Schiefe der Ekliptik 
  T=MJD+2400000.5;
  [xRA,xDec]=KeplerSonne(T,eps);  % Keplerlösung f. Rektaszension/Deklination 
  tau=GMST(MJD)+deg2rad(lambda)-xRA;
  tau = wrapToPi(tau);          % Berechnung des Stundenwinkels 
  sinusalt=sinphi*sin(xDec)+cosphi*cos(xDec)*cos(tau);
end

