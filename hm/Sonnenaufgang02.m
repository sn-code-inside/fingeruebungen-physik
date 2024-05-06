% -------------------------------------------------------------------------
% Sonnenaufgang02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die Aufgangszeiten auf Basis der 
% Keplergleichung 
% a) für die Erdbahn mit einer iterativen Methode (Quadratische 
% Interpolation) zur Suche der Nullstellen für Auf- und
% Untergänge und 
% b) zum Vergleich auf Basis der ZGL-Näherung.
% Die Abweichungen werden dargestellt.
%
% Achten Sie bei der Ausführung des Programms auf Ihre Kommandozeile
% (Command Window) und bestätigen Sie mit "J", falls die Daten für
% Oberkochen verwendet werden sollen. Möchten Sie dieses Programm für einen
% anderen Ort verwenden, wählen Sie "N".
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
pi2 = 2*pi;                   % 2*pi
eps0=deg2rad(23.43929111);
hr=deg2rad(-50/60);

% Abfrage der Eingabedaten in der Kommandozeile
Oko     = input('\Oberkochen J/N  :                 ... ','s');
if Oko == 'J' || Oko =='j' 
    phiOko=48.7854368;
    lambdaOko=10.0990193;
    Zone=1;
    phi=phiOko;
    lambda=lambdaOko;
else
    phi     = input('\Geographische Breite +-BB.bb      ... ');
    lambda  = input('\Geographische Laenge +-LL.ll      ... ');
    Zone    = input('\Zeitzone inkl. DST +-UT           ... ');
end

phistr=sprintf(' %+4.1f°',phi);
lstr=sprintf(' %+4.1f°',lambda);    

phi=deg2rad(phi);
sinphi=sin(phi);
cosphi=cos(phi);

%_________________________________________________________________________

% Hauptprogramm
%
%------------------------------------------------------------------------------

MOZ= lambda/15;  % mittlere Ortszeit (s. Kap. Raum und Zeit)
eps0=deg2rad(23.43929111);      % Schiefe der Ekliptik

dt1 = datetime('2000-01-01 12:00:00');
T = juliandate(dt1) - Zone/24; % Julianisches Datum
MJuDa = juliandate(dt1,'modifiedjuliandate')- Zone/24;% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');


T1 = T-1;
MJuDa1=MJuDa-1;
MJuDa0=MJuDa-1.5;

%Felder für die Näherung
sunrise0 = zeros(1,366);
sunset0  = zeros(1,366);
zeitdiff0= zeros(1,366);
kulmi0   = zeros(1,366);
tauw0    = zeros(1,366);

%Felder für die numerische Lösung
sunrise1 = zeros(1,366);
sunset1  = zeros(1,366);
zeitdiff1= zeros(1,366);
kulmi1  = zeros(1,366);
tauw1    = zeros(1,366);

%-----------------------------------------------------------------------------

% Schleife ueber 366 aufeinanderfolgende Tage
  
T_vector=(T+1):(T+366);
MJuDa_V1 = (MJuDa1+1):(MJuDa1+366);
MJuDa_V0 = (MJuDa1+0):(MJuDa0+366);
% ZGL Näherung 
%Berechnung nach Formel mit Interpolation der Deklination und RA
[RA_vector,Dec_vector]=KeplerSonne(T_vector,eps0); 
tauw0=(GMST(MJuDa_V1)+deg2rad(lambda)-RA_vector); 
tauw0 = rad2deg(tauw0)/15 -MOZ;
tauw0 = 24*wrapToPi(pi2*tauw0/24)/pi2;
q=(sin(hr)-sin(phi)*sin(Dec_vector))./(cos(phi)*cos(Dec_vector));
kulmi0=12-tauw0-MOZ;  %Kulmination

% Abfangen von zirkumpolaren und sub-Horizont Situationen
for k=1:366
    if abs(q(k)) < 1
        zeitdiff0(k)   = acosd(q(k))/15; 
        sunrise0(k)    = kulmi0(k) - zeitdiff0(k);
        sunset0(k)     = kulmi0(k) + zeitdiff0(k);
        taglaenge(k)   = 2*zeitdiff0(k);
    else
        if q(k)<0 
            taglaenge0(k)= 23.99;
         else
            taglaenge0(k)= 0.05; 
            kulmi0(k)=NaN;  %Kulmination
        end
        zeitdiff0(k) = acosd(q(k))/15; 
        sunrise0(k) = NaN; 
        sunset0(k)  = NaN; 
    end
end   
daylength1= sunset0-sunrise0;

% Schleife ueber 366 aufeinanderfolgende Tage
% Numerische Berechnung mit quadratischer Interpolation
kulmi1 = kulmi0;
for k=1:366
 T=T+1;
 MJuDa1=MJuDa1+1;
 MJuDa0=MJuDa0+1;
 %Berechnung mit quadratischer Interpolation
 [sunrisex, sunsetx, xRise, xSet, xAbove] = FindRiseSet(MJuDa0,lambda,phi);
 if xRise
    sunrise1(1,k) = (wrapTo360((sunrisex)*15))/15;
 else
    sunrise1(1,k) = NaN;
 end
 if xSet
    sunset1(1,k) = (wrapTo360((sunsetx)*15))/15;
 else
    sunset1(1,k) = NaN;
 end
 if xAbove==-1
    kulmi1(k) = NaN;
 end
end
daylength1= sunset1-sunrise1;

%------------------------------------------------------------------------------
% %  Ausgabe

out1=sunrise1;
out2=sunset1;
out3=kulmi1;
out4=daylength1;
out5(:) = sunrise1(1,:)- sunrise0(1,:);
out6(:) = sunset1(1,:) - sunset0(1,:);

% 'SA','SU','Mittag','Taglänge';
figure();
for k=1:12
    XDataTick(k) = datetime(2020,k,1);
end
XDataTick(13) = datetime(2021,1,1);
xData=datetime('2019-12-31') + caldays(1:366);
plot(xData, out1,'Color',Colors(2,:),'LineWidth',2); 
hold on;
plot(xData, out2,'Color',Colors(3,:),'LineWidth',2); 
plot(xData, out3,'Color',Colors(4,:),'LineWidth',2); 
plot(xData, out4,'Color',Colors(5,:),'LineWidth',2); 
grid on;
ax = gca;
ax.XTick = XDataTick;
datetick('x','mmm','keepticks');
lgd=xlabel('Tag');
lgd.FontSize=12;
lgd=legend('SA','SU','Mittag','Taglänge','location','northeast','NumColumns',2);
lgd.FontSize=16;
legend boxoff;
txt = ['SA,SU,Taglänge'];
text(5,22.66, txt,'FontSize',16);
txt= ['Breite {\phi} = ',phistr,' und Länge {\lambda} = ',lstr];
text(5,21.5,txt,'FontSize',16);
ylim([0,24]);
lgd=ylabel('Zeit in MOZ bzw. Länge in h');
lgd.FontSize=12;
yticks manual;
yticks([0 3 6 9 12 15 18 21 24]);
yticklabels({'00:00' '03:00' '06:00' '09:00' '12:00' '15:00' '18:00' '21:00' '24:00'});
hold on;
set(gca,'FontSize',18);

%Fehleranzeige nur für Bereiche ohne Polarsommer/Polarwinter
if abs(rad2deg(phi)) < (90 - rad2deg(eps0) - 0.1)  
    figure();
    for k=1:12
        XDataTick(k) = datetime(2020,k,1);
    end
    XDataTick(13) = datetime(2021,1,1);
    xData=datetime('2020-01-01') + caldays(1:366);
    plot(xData,out5*3600,'Color',Colors(2,:),'LineWidth',2);
    hold on;
    plot(xData,out6*3600,'Color',Colors(3,:),'LineWidth',2);
    ax = gca;
    ax.XTick = XDataTick;
    datetick('x','dd.mm.','keepticks');
    lgd=title('Abweichung numerische Lösung von ZGL-Näherung');
    lgd.FontSize=18;
    grid on;
    lgd=legend('SA','SU','location','northeast');
    lgd.FontSize=18;
    legend boxoff;
    lgd=xlabel('Tag');
    lgd.FontSize=18;
    lgd=ylabel('Abweichung in sec');
    lgd.FontSize=18;
    set(gca,'FontSize',18);
end

fprintf('\n Alle Zeiten in MOZ \n');
for n=1:366
     MJuDa1=MJuDa0+n+1;
     DateBer = string(datetime(MJuDa1,'ConvertFrom','modifiedjuliandate'),'dd.MM.yyyy');
     SunRiseStr  = StrHMS(sunrise1(n));
     SunSetStr   = StrHMS(sunset1(n));
     KulmiStr    = StrHMS(kulmi1(n));
     if isnan(daylength1(n)) 
         if isnan(kulmi1(n))
            DayLength = '00:00:00';
         else
            DayLength = '24:00:00';
         end
     else
        DayLength   = StrHMS(daylength1(n));
     end    
     fprintf('Tag: %s | SA: %s | SU: %s | K: %s | Taglänge: % s \n', DateBer, SunRiseStr, SunSetStr, KulmiStr, DayLength);
end

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
