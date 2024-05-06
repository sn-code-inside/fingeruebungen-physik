% -------------------------------------------------------------------------
% VenusTransit02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Koordinaten der Venus und des Merkur bei Transits.
% Die Zeiten muessen naeherungsweise z.B. aus InnerePlanetenPQR.m bekannt
% sein. Daraus wird die Astronomische Einheit aus den Daten eines 
% Venustransits bestimmt.
% a) Parallaxen-Methode: Die beobachteten Werten fuer zwei geographischen
% Orten zu einer festen Zeit(UT) werden dem Internet entnommen, daraus 
% wird die Parallaxe bestimmt und daraus dann die AE berechnet.  
% b) Kontaktzeitmethode durch Beobachtung an einem Standort.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

%%
% -------------------------------------------------------------------------
%  a) PARALLAXENMETHODE nach Halley 
%

% Koordinaten der beiden Orte Tromsoe und Canberra
RE = 6378;
phi(1)= 35.28 ;
lambda(1)= 18.95;
phi(2) = -35.28;
lambda(2) = 149.13;
dt = datetime('2012-06-06 03:12:15');
T1 = juliandate(dt);
bw  = deg2rad(0.5317);    % Winkeldurchmesser der Sonne am 8.6.2004
dnu = 0.021*bw; % Winkelabstand der Transitpunkte
gamma = 1.382;  % Verhaeltnis Erd-Venus-Abstand
MJuDa = juliandate(dt,'modifiedjuliandate');
% Sternzeiten der beiden Orte
theta(1) = GMST(MJuDa) + deg2rad(lambda(1));
theta(2) = GMST(MJuDa) + deg2rad(lambda(2));
% Koordinaten der Sonne
epsErad=deg2rad(EpsErde(T1));
[alpha,delta]=KeplerSonne(T1, epsErad);
r_Sun = CalcXYZfromAngles([1; alpha; delta]);
% Koordinaten der Orte
rT=CalcXYZfromAngles([1; theta(1); deg2rad(phi(1))]);
rC=CalcXYZfromAngles([1; theta(2); deg2rad(phi(2))]);
verb=rT-rC;
y=sum(verb.*r_Sun)/norm(verb)/norm(r_Sun);
w=rad2deg(acos(y));
d_eff= RE*sqrt(sum(verb.*verb))*sind(w);
AEge = d_eff/(gamma-1)/dnu;
PI = RE/AEge;
fprintf('\n');
fprintf('Effektiver Abstand von \n')
fprintf('Tromsoe und Canberra : %12.1f km \n \n', d_eff);
fprintf('Astronomische Einheit: %12.1f km \t ', AEge);
fprintf('Parallaxe: %s \n', StrDMS(rad2deg(PI)));


%%
% -------------------------------------------------------------------------
%  b) KONTAKTZEITMETHODE 
%
% Geozentrische Koordinaten von Oberkochen zu Beginn am Ende des Transits 
RE = 6378;
phi(1)= 48.78 ;
lambda(1)= 10.1;
dt1 = datetime('2004-06-08 05:21:30');
T1 = juliandate(dt1);
dt2 = datetime('2004-06-08 11:15:00');
T2 = juliandate(dt2);
% Sternzeiten der beiden Kontakte
MJuDa1 = juliandate(dt1,'modifiedjuliandate');
MJuDa2 = juliandate(dt2,'modifiedjuliandate');
% Sternzeiten der beiden Orte
theta(1) = GMST(MJuDa1) + deg2rad(lambda(1));
theta(2) = GMST(MJuDa2) + deg2rad(lambda(1));
% Koordinaten der Sonne
eps=deg2rad(23.43929111);
[alpha,delta]=KeplerSonne((T1+T2)/2, eps);
ekl_Sun = CalcXYZfromAngles([1; alpha; delta]);
ekl1=CalcXYZfromAngles([1; theta(1); deg2rad(phi(1))]);
ekl2=CalcXYZfromAngles([1; theta(2); deg2rad(phi(1))]);
verb=ekl1-ekl2;
y=sum(verb.*ekl_Sun)/norm(verb)/norm(ekl_Sun);
w=rad2deg(acos(y));
d_eff= RE*sqrt(sum(verb.*verb))*sind(w);

bw = 0.00928; % Winkeldurchmesser der Sonnne von der Erde aus
TV = 224.701; % siderische Umlaufzeit Venus
TE = 365.256; % siderische Umlaufzeit Erde
gamma  = (TE^2/TV^2)^(1/3);
% Basis Winkelgeschwindigkeit der Projektion
OmegaB = 2*pi*(1/TV -1/TE)/86400/(gamma-1); 
T_max  = bw/OmegaB;
ls     = 0.775;
T_1    = ls*T_max;
T_TR   = 21780;
DelT   = linspace(-120,120,240);
DA_A   = OmegaB.*DelT*100./(bw*ls-OmegaB*(T_TR+abs(DelT)));

d_effv = [d_eff, 4000, 5000, d_eff, d_eff];
lsv    = [ls, ls, ls, 0.995*ls, 1.005*ls];

T = linspace(21000,22000,1000);
AE(1,:)  = d_effv(1)/(gamma-1)./(bw*ls-OmegaB*T); 
AE(2,:)  = AE(1,:)*d_effv(2)/d_eff;
AE(3,:)  = AE(1,:)*d_effv(3)/d_eff;
AE(4,:)  = d_eff/(gamma-1)./(bw*lsv(4)-OmegaB*T);
AE(5,:)  = d_eff/(gamma-1)./(bw*lsv(5)-OmegaB*T);

for k=1:5 
   lgdtitle(k,:) = strjoin(["\it d\rm_{eff} = ", ...
       num2str(round(d_effv(k)),'%4d'),"\it l\rm_S = ", ...
       num2str(lsv(k),'%4.3f')]);
end

% -------------------------------------------------------------------------
%%

% Graphische Ausgabe

lgdtitle(6,:) = ' Astron. Einheit';
lgdtitle(7,:) = ' Messwerte';
Messwert0 = 21780;
Messwertminus = 21720;
Messwertplus = 21840;

figure('Name', "AE Bestimmung ueber Kontaktzeitmethode"),
subplot(1,2,1);
plot(DelT,DA_A,'Color', Colors(2,:),'Linewidth',2);
xlim([-120 120]);
ylim([-50 50]);
xlabel('Zeitfehler in s ')
ylabel('rel. Messfehler in %');
grid on;
set(gca,'FontSize',16);

AE0 = 149.6 + 0*T;
subplot(1,2,2);
plot(T/1000,AE(1,:)/1e6,'Color', Colors(2,:),'Linewidth',2);
hold on;
plot(T/1000,AE(2,:)/1e6,'Color', Colors(3,:),'Linewidth',1);
plot(T/1000,AE(3,:)/1e6,'Color', Colors(4,:),'Linewidth',1);
plot(T/1000,AE(4,:)/1e6,'Color', Colors(3,:),'Linewidth',1,'LineStyle',...
    '-.');
plot(T/1000,AE(5,:)/1e6,'Color', Colors(4,:),'Linewidth',1,'LineStyle',...
    '-.');
plot(T/1000,AE0,'Color', Colors(3,:),'Linewidth',1);
y = linspace(50,250,100);
line(Messwert0/1000+y*0,y,'Color','k','Linewidth',2,'LineStyle',':');
line(Messwertplus/1000+y*0,y,'Color','k','LineStyle',':','Linewidth',2);
line(Messwertminus/1000+y*0,y,'Color','k','LineStyle',':','Linewidth',2);
ylim([50 200]);
xlim([21.1 21.9]);
legend(lgdtitle,'location','northwest');
legend box off;
grid on;
grid minor;
xlabel('\it T\rm_{TR} in 1000 s ')
ylabel('AE in Mio. km');
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
