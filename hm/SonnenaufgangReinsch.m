% -------------------------------------------------------------------------
% SonnenaufgangReinsch.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Analytische Formel für Sonnenaufgang und Vergleich mit numerischer Berechnung
% Formel nach M. Reinsch, Am. J. Phys. 71, 1242 (2003)
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
ups0 = 0;

pi2=2*pi;
eps=deg2rad(23.43929111);
phiOko=48.7854368;
phi=phiOko;
phi=deg2rad(phi);
phistr=sprintf(' %+4.1f°',rad2deg(phi));

T_sun = 86400;
T_star= 86184;
Delta0 = pi2/T_star - pi2/T_sun;

t=linspace(1, 366, 366);
% Zeit-Parameter
dt1 = datetime('2000-01-01 12:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');

% Frühlingspunkt
tFP = datetime(2000,03,21);
FP=day(tFP,'dayofyear');
T0=Jd(MJuDa1)-1;
MJuDa0=MJuDa1-1.5;
% Frühlingspunkt
T=(T0):(T0+366);
MJuDa = (MJuDa1+FP):(MJuDa1+366+FP);
T=(T-2451545.0)/36525;

% Berechnung über analytische Formel 
% Formel nach M. Reinsch, Am. J. Phys. 71, 1242 (2003)
% Mittlere Anomalie und ekliptikale Laenge
M  = pi2 * Frac(0.993133*0  + 99.997361*T); 
L  = pi2 * Frac(0.7859453*0 + M/pi2+(6893*sin(M)+72*sin(2*M)+0.0*sin(3*M)+0*T)/(1296.0*1000)) ;
Delta = GMST(MJuDa);
X  = sqrt((cos(eps)*sin(L).*sin(Delta)+cos(L).*cos(Delta)).^2 + (cos(eps)*sin(L).*cos(Delta)-cos(L).*sin(Delta)).^2);
xi = -atan2(cos(eps)*sin(L).*cos(Delta)-cos(L).*sin(Delta),(cos(eps)*sin(L).*sin(Delta)+cos(L).*cos(Delta)));
T_SA =12-(acos((-sin(L)*sin(eps)*tan(phi))./X)-xi)*12/pi;
T_SU =12+(acos((-sin(L)*sin(eps)*tan(phi))./X)-xi)*12/pi;

% Shift auf Jahresanfang
for k=1:(366-FP)
  H_SA(k+FP)=T_SA(k);
  H_SU(k+FP)=T_SU(k) ; 
end
for k=1:FP
  H_SA(k)=T_SA(366-FP+k);
  H_SU(k)=T_SU(366-FP+k); 
end

% Numerische Berechnung mit quadratischer Interpolation
for k=1:366
 T=T+1;
 MJuDa1=MJuDa1+1;
 MJuDa0=MJuDa0+1;
% Quadratischer Interpolation
 [sunrisex, sunsetx, xRise, xSet, xAbove] = FindRiseSet(MJuDa0,0,phi);
 if xRise
    sunrise0(1,k) = (wrapTo360((sunrisex)*15))/15;
 else
    sunrise0(1,k) = NaN;
 end
 if xSet
    sunset0(1,k) = (wrapTo360((sunsetx)*15))/15;
 else
    sunset0(1,k) = NaN;
 end
 if xAbove<0
    kulmi1(k) = NaN;
 end
end

%------------------------------------------------------------------------------
%%  Ausgabe

out1=sunrise0(1,:);
out2=sunset0(1,:);

figure();
xData=datetime('2019-12-31') + caldays(1:366);
plot(xData, out1,'Color',Colors(2,:),'LineWidth',1); 
hold on;
plot(xData, out2,'Color',Colors(3,:),'LineWidth',1); 
plot(xData, H_SA, 'Color', Colors(2,:),'LineStyle','-.','LineWidth',2);
plot(xData, H_SU, 'Color', Colors(3,:),'LineStyle','-.','LineWidth',2);
grid on;
ax = gca;
for k=1:12
    XDataTick(k) = datetime(2020,k,1);
end
XDataTick(13) = datetime(2021,1,1);
ax.XTick = XDataTick;
datetick('x','mmm','keepticks');
lgd=xlabel('Tag');
lgd.FontSize=12;
lgd=legend('SA analyt.','SU analyt.','SA numer.','SA numer.','location','northeast','NumColumns',2);
lgd.FontSize=16;
legend boxoff;
txt= ['Breite {\phi} = ',phistr];
text(5,21.5,txt,'FontSize',16);
ylim([0,24]);
lgd=ylabel('Zeit in UT');
lgd.FontSize=12;
yticks manual;
yticks([0 3 6 9 12 15 18 21 24]);
yticklabels({'00:00' '03:00' '06:00' '09:00' '12:00' '15:00' '18:00' '21:00' '24:00'});
hold on;
set(gca,'FontSize',14);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------