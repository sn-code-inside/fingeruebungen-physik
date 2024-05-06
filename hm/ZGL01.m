% -------------------------------------------------------------------------
% ZGL01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die ZGL auf Basis einer (Mittelpunktsgleichung 
% basierten) Keplerlösung für die Sonnenbahn und vergleicht diese mit
% empirischen Werten und den Berechnungen des JPL.
% (Äquinoktium und Ekliptik des Datum J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
dt1 = datetime('2000-01-01 12:00:00');
T = juliandate(dt1);
MJuDa = juliandate(dt1,'modifiedjuliandate');
DatumStr1 =  string(dt1,'dd.MM.yyyy');

T0=(T-2451545)/36525;
eps=deg2rad(23.43929111-(46.8150+0.00059*T0-0.001813*T0*T0)*T0/3600);
pi2 = 2*pi;

dt2 = datetime('2000-12-31 12:00:00');
DatumStr2 =  string(dt2,'dd.MM.yyyy');


%-------------------------------------------------------------------------
%Beginn Rechnung

T_vector=T:(T+365);
MJuDa_vector = MJuDa:(MJuDa+365);
Timey = 1:366;
Null = zeros(1,366);
% ZGL über Kepler-Lösung (Mittelpunktsgleichung)
% Berechnung nach Keplerlösung
[RA_vector,Dec_vector]=KeplerSonne(T_vector,eps); 
% Berechnung des Stundenwinkels
tau0 = rad2deg((GMST(MJuDa_vector)-RA_vector))/15;
%Stundenwinkel in Zeitminuten
tau0 = 24*60*wrapToPi(pi2*tau0/24)/pi2; 

% Empirische Daten Klausmann (100 Jahre gemittelt)
fid=fopen('ZGLempirisch.dat','r');
if fid ==1 
    disp('File open not successful')
else
   Emp_ZGL =-fscanf(fid,'%f', [1, 366]);
end
closeresult =fclose(fid);
if closeresult ==0
%     disp('Data successfully loaded');
else
     disp('Dataload not successful');
end
Delta1= (tau0-Emp_ZGL)*60;

% ZGL über Horizons JPL 
fid=fopen('ZGLhorizonJPL.dat','r');
if fid ==1 
    disp('File open not successful')
else
   Emp_ZGL_Alpha =fscanf(fid,'%f', [1, 366]);
end
closeresult =fclose(fid);
if closeresult ==0
%     disp('Data successfully loaded');
else
    disp('Dataload not successful');
end

% Berechnung des Stundenwinkels
Emp_ZGL_Alpha =deg2rad(Emp_ZGL_Alpha);
tau1 = rad2deg((GMST(MJuDa_vector)-Emp_ZGL_Alpha))/15;
% Stundenwinkel in Zeitminuten
Emp_ZGL_Alpha = 24*60*wrapToPi(pi2*tau1/24)/pi2; 
Delta2= (tau0-Emp_ZGL_Alpha)*60;
for i=2:366
   if abs(Delta2(i)) > 10
      Delta2(i)=Delta2(i-1);   
   end
end


%------------------------------------------------------------------------------
% Graphische Ausgabe

header1='Zeitgleichung'; 
figure('Name',header1);
subplot(2,1,1);
plot(Timey,tau0,'Color',Colors(3,:),'LineWidth',2);
hold on;
plot(Timey,Null,'Color','k');
ylim([-20 20]);
xlim([0 366]);
grid on;
header2=sprintf('%s %s %s %s','Epoche :', DatumStr1,' - ',DatumStr2);
title(header2);
xlabel('Tag')
ylabel('ZGL Minuten');
legend(header2,'location','northwest','NumColumns',2);
legend boxoff;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

subplot(2,1,2);
plot(Timey,Delta1,'Color', Colors(3,:),'LineWidth',2);
hold on;
plot(Timey,Delta2,'Color', Colors(2,:),'LineWidth',2);
plot(Timey,Null,'k');
ylim([-20 20])
xlim([1 366]);
header2 = strcat( 'Abweichung von empirischen Daten');
% title(header2);
xlabel('Tag')
ylabel('Abweichung Sekunden');
ax=gca;
ax.FontSize=18;
legend(' Klausmann',' JPL Horizon','location','northwest','NumColumns',2);
legend boxoff;
grid on;
% grid minor;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
