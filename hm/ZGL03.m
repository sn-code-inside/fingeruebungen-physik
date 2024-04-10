% -------------------------------------------------------------------------
% ZGL03.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet die ZGL auf Basis einer (auf der Mittelpunktsgleichung 
% basierten) Keplerlösung für die Erdbahn und vergleicht diese mit 
% verschiedenen publizierten Näherungsformeln und den Berechnungen des JPL.  
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
DatumStr =  string(dt1,'dd.MM.yyyy');
T0=(T-2451545)/36525;
eps=deg2rad(23.43929111-(46.8150+0.00059*T0-0.001813*T0*T0)*T0/3600);
pi2 = 2*pi;
U = 365.25;
DeltaR = 5;
RAE = 149.6;

%-------------------------------------------------------------------------
%Begin Rechnung

T_vector=(T):(T+365);
MJuDa_vector = (MJuDa):(MJuDa+365);
Timey = 1:366;

% ZGL über Kepler-Lösung (Mittelpunktsgleichung) numerisch
[RA_vector,Dec_vector]=KeplerSonne(T_vector,eps); %Berechnung nach Keplerlösung
tau0 = rad2deg((GMST(MJuDa_vector)-RA_vector))/15;
tau0 = 24*60*wrapToPi(pi2*tau0/24)/pi2; %Stundenwinkel = ZGL0 in Zeitminuten
tauw(1,:)=tau0;

% ZGL nach Feynmann
% tauw(2,:) = -9.863*sin(0.03443*(Timey)+0.344284)-7.6599*sin(0.017202*(Timey)-0.051643);

% ZGL nach Formel (3.83)
Timey1=Timey;
tauw(2,:) = -7.364186*sin(0.017202*(Timey1)-0.062882);
tauw(2,:) = tauw(2,:)-9.93481*sin(0.034404*(Timey1)+0.361835);
tauw(2,:) = tauw(2,:)-0.329614*sin(0.051606*(Timey1)+0.321940);
tauw(2,:) = tauw(2,:)-0.21222*sin(0.068808*(Timey1)+0.730548);

tauw(2,:) = -7.3635*sin(0.017202*(Timey1)-0.06266);
tauw(2,:) = tauw(2,:)-9.9348*sin(0.034404*(Timey1)+0.361835);
tauw(2,:) = tauw(2,:)-0.32961*sin(0.051606*(Timey1)+0.322223);
tauw(2,:) = tauw(2,:)-0.21222*sin(0.068808*(Timey1)+0.730697);

% ZGL nach Astro-Lexikon bzw. Wetzel2007 (Internet www.astronomie.info)
tauw(3,:) = -10.26*sin(0.0337*(Timey)+0.465)-7.794*sin(0.01787*(Timey)-0.168);

% ZGL nach Rempel 11/2019 ab initio
taun     = pi2*Modulo((Timey+10)/U,1);
Deltatau = (DeltaR/RAE)*sin(pi2*(Timey-3)/U);
omegan      = taun+Deltatau;
Arctan = atan2(sin(omegan),(cos(eps)*cos(omegan)));
Deltaomega  = mod(Arctan+pi2,pi2);
tauw(4,:) = -(Deltaomega - taun)*24*60/pi2/(1+1/U);

fid=fopen('ZGLempirisch.dat','r');
if fid ==1 
    disp('File open not successful')
else
    Emp_ZGL =-fscanf(fid,'%f ', [1, inf]);
end
closeresult =fclose(fid);
if closeresult ==0
%     disp('Data successfully loaded');
else
    disp('Dataload not successful');
end

% ZGL über Horizons JPL (=exakte Lösung)
fid=fopen('ZGLhorizonJPL.dat','r');
if fid ==1 
    disp('File open not successful')
else
   JPL_ZGL_Alpha =fscanf(fid,'%f', [1, inf]);
end
closeresult =fclose(fid);
if closeresult ==0
%     disp('Data successfully loaded');
else
    disp('Dataload not successful');
end
% Berechnung des Stundenwinkels aus JPL Daten
JPL_ZGL_Alpha =deg2rad(JPL_ZGL_Alpha);
tau1 = rad2deg((GMST(MJuDa_vector)-JPL_ZGL_Alpha))/15;
%Stundenwinkel in Zeitminuten
JPL_ZGL_Alpha = 24*60*wrapToPi(pi2*tau1/24)/pi2; 
Delta= (tau0-JPL_ZGL_Alpha)*60;
tauw(5,:)=JPL_ZGL_Alpha;
for i=2:365
   if abs(Delta(i)) > 10
      Delta(i)=Delta(i-1);  
      tauw(5,i)=0.5*(JPL_ZGL_Alpha(i-1)+JPL_ZGL_Alpha(i+1));
   end
end

%------------------------------------------------------------------------------
% Graphische Ausgabe


titlestr   =  strings([5,25]);
titlestr(1,:)= ' Numerisch (MPG)'; %einfache Keplerlösung
titlestr(2,:)= ' Formel (3.83) ';
titlestr(3,:)= ' Rempel2019';
titlestr(4,:)= ' Wetzel2007';
titlestr(5,:)= ' JPL';  %JP = exakt

header1='Zeitgleichung aus verschiedenen Quellen ';
figure('Name',header1);
for k=1:5
    plot(Timey, tauw(k,:), 'Color', Colors(k,:));
    hold on;
end
ylim([-30 30]);
xlim([0 370]);
grid on;
grid minor;
title(header1);
xlabel('Tag')
ylabel('ZGL in min');
legend(titlestr(1),titlestr(2),titlestr(3),titlestr(4),titlestr(5),'location','southeast');
legend boxoff;


TimeStr=string(hours(12)+minutes(0)+seconds(0),'hh:mm');
header1='Abweichung Zeitgleichungen von Empirischen Daten'; 
figure('Name',header1);
for k=1:5
    plot(Timey, (tauw(k,:)-Emp_ZGL)*60, 'Color', Colors(k,:));
    hold on;
end
ylim([-30 30]);
xlim([1 370]);
grid on;
grid minor;
title(header1);
xlabel('Tag')
ylabel('Abweichung in s');
legend(titlestr(1),titlestr(2),titlestr(3),titlestr(4),titlestr(5),'location','south');
legend boxoff;

header1='Abweichung verschiedener Zeitgleichungen von der VSOP Lösung';
figure('Name',header1);
for k=1:4
    subplot(2,2,k);
    plot(Timey, (tauw(k,:)-tauw(5,:))*60, 'Color', Colors(k,:));
    hold on;
    title(titlestr(k));
    ylim([-20 20]);
    xlim([1 370]);
    xlabel('Tag')
    ylabel('Abweichung in s');
    grid on;
    grid minor;
end
header1='Abweichung verschiedener Zeitgleichungen von der JPL-VSOP Lösung';
figure('Name',header1);
for k=1:4
    plot(Timey, (tauw(k,:)-tauw(5,:))*60, 'Color', Colors(k,:),'LineWidth',2);
    hold on;
    header2=titlestr(k); 
end
ylim([-60 60]);
xlim([1 370]);
title(header1);
xlabel('Tag')
ylabel('Abweichung in s');
grid on;
grid minor;
legend(titlestr(1),titlestr(2),titlestr(3),titlestr(4),'location','south');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
