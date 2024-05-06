% -------------------------------------------------------------------------
% Analemma01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Annalemma über die ZGL die auf einer auf der 
% Mittelpunktsgleichung basierten Keplerlösung für die Sonnenbahn beruht. 
% Die Exzentrizität und die Ekliptikschiefe werden variiert. 
% Es wird die ZGL für UT=12:00 berechnet. 
% (Äquinoktium und Ekliptik J2000) 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
titlestr   =  strings([6,1]);

dt = datetime('2000-01-01 12:00:00');
T = juliandate(dt);
MJuDa = juliandate(dt,'modifiedjuliandate');
DatumStr =  string(dt,'dd.MM.yyyy');

%für Schaltjahre (2000 war ein solches)
BeginnMonat=[1 32 61 92 122 153 183 214 245 275 306 336];
%für Nicht-Schaltjahre
%BeginnMonat=[1 32 60 91 121 152 182 213 244 274 305 335];

% T=Jd(MJuDa);                %Julianisches Datum
T0=(T-2451545)/36525;         %Zeit seit 1.1.2000 12:00

%zeitabhängige Berechnung der Ekliptikschiefe
eps0=23.43929111-(46.8150+0.00059*T0-0.001813*T0*T0)*T0/3600;

e0 =0.016709;               %Exzentrizität Bsp.: Erde 0.016709
pi2 = 2*pi;

%-------------------------------------------------------------------------
%Beginn Rechnung
T_vector= T:(T+365);
MJuDa_vector= MJuDa:(MJuDa+365);
Timey=1:366;
Null = zeros(1,366);

% ZGL über Kepler-Lösung (Mittelpunktsgleichung) Variation eps
for m=1:3
    eps=deg2rad(eps0+(m-2)*5);
    %Berechnung nach Keplerlösung
    [RA_vector,Dec_vector]=KeplerSonne(T_vector,eps);
    % Berechnung des Stundenwinkels im Zeitmaß (Zeitstunden)
    tau0 = rad2deg((GMST(MJuDa_vector)-RA_vector))/15;
    %Stundenwinkel in Zeitminuten symmetrisch zu Null
    tau0 = 24*60*wrapToPi(pi2*tau0/24)/pi2;
    tauw(m,:)=tau0;
    Deklination(m,:)= rad2deg(Dec_vector);
    titlestr(m,:)= sprintf('eps = %2.1f °',rad2deg(eps));
end

% ZGL über Kepler-Lösung (Mittelpunktsgleichung) Variation Exzentrizizät
% exz ist der Faktor, mit dem wir die Exzentrizität e0 multiplizieren
exz=[0.0001,1,2];
for m=1:3
    eps=deg2rad(eps0);
    exzx= exz(m);
    %Berechnung nach Keplerlösung
    [RA_vector,Dec_vector]=KeplerSonneEx(T_vector,eps,exzx);
    % Berechnung des Stundenwinkels im Zeitmaß (Zeitstunden)
    tau0 = rad2deg((GMST(MJuDa_vector)-RA_vector))/15;
    %Stundenwinkel in Zeitminuten symmetrisch zu Null
    tau0 = 24*60*wrapToPi(pi2*tau0/24)/pi2;
    tauw(m+3,:)=tau0;
    Deklination(m+3,:)= rad2deg(Dec_vector);
    titlestr(m+3,:)= sprintf('e = %5.4f',exzx*e0);
end

%------------------------------------------------------------------------------
% Graphische Ausgabe
header1='Zeitgleichung'; 
figure('Name',header1);
subplot(1,2,1);
plot(Timey, tauw(1,:),'Color',Colors(2,:),'LineWidth',2);
hold on;
plot(Timey, tauw(2,:),'Color',Colors(3,:),'LineWidth',2);
plot(Timey, tauw(3,:),'Color',Colors(4,:),'LineWidth',2);
ylim([-30 30]);
xlim([0 370]);
header2=strcat('Variation der Ekliptikschiefe');
title(header2);
grid on;
xlabel('Tag')
ylabel('ZGL in min');
legend(titlestr(1),titlestr(2),titlestr(3),'location','best');

subplot(1,2,2);
plot(Timey, tauw(4,:),'Color',Colors(2,:),'LineWidth',2);
hold on;
plot(Timey, tauw(5,:),'Color',Colors(3,:),'LineWidth',2);
plot(Timey, tauw(6,:),'Color',Colors(4,:),'LineWidth',2);
ylim([-30 30]);
xlim([0 370]);
header2=strcat('Variation der Exzentrizität');
title(header2);
grid on;
xlabel('Tag')
ylabel('ZGL in min');
legend(titlestr(4),titlestr(5),titlestr(6),'location','best');

header1='Analemma'; 
figure('Name',header1);
subplot(1,2,1);
plot(tauw(1,:),Deklination(1,:),'Color',Colors(2,:),'LineWidth',2);
hold on;
plot(tauw(2,:),Deklination(2,:),'-+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(3,:),'LineStyle','-');
plot(tauw(3,:),Deklination(3,:),'Color',Colors(4,:),'LineWidth',2);
for k=1:12
        mylabels(k,:)=sprintf('1.%02u',k);
        xL(k)=tauw(2,BeginnMonat(k));
        yL(k)=Deklination(2,BeginnMonat(k));
        LabelPoints(xL, yL ,mylabels,'w',0.2,0,'FontSize',12,'Color',Colors(3,:));
end

ylim([-40 30]);
xlim([-25 25]);
grid on;
header2=strcat('Variation der Ekliptikschiefe');
% title(header2);
xlabel('Zeit in min')
ylabel('Deklination °');
legend(titlestr(1),titlestr(2),titlestr(3),'location','southeast');
legend boxoff;
set(gca,'FontSize',16);

subplot(1,2,2);
plot(tauw(4,:),Deklination(4,:),'Color',Colors(2,:),'LineWidth',2);
hold on;
plot(tauw(5,:),Deklination(5,:),'-+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(3,:),'LineStyle','-');
plot(tauw(6,:),Deklination(6,:),'Color',Colors(4,:),'LineWidth',2);
for k=1:12
        mylabels(k,:)=sprintf('1.%02u',k);
        xL(k)=tauw(5,BeginnMonat(k));
        yL(k)=Deklination(5,BeginnMonat(k));
        LabelPoints(xL, yL ,mylabels,'w',0.2,0,'FontSize',12,'Color',Colors(3,:));
end
ylim([-40 30]);
xlim([-25 25]);
grid on;
header2=strcat('Variation der Exzentrizität');
% title(header2);
xlabel('Zeit in min ')
ylabel('Deklination °');
legend(titlestr(4),titlestr(5),titlestr(6),'location','southeast');
legend boxoff;
    
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
