% -------------------------------------------------------------------------
% Besselreihe.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Exzentrische Anomalie auf Basis der Besselreihe
% Berechnung von Bahnkurven nach Besselreihe
% Betrachtun zur Genauigkeit der Besselreihe
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
dt1 = datetime('2000-01-01 00:00:00');
T1 = juliandate(dt1);
dt2 = datetime('2001-01-01 00:00:00');
T2 = juliandate(dt2);
MJuDa = juliandate(dt1,'modifiedjuliandate');
DatumStr1 =  string(dt1,'dd.MM.yyyy');

% Annahme: 1 Jahr = 365 Tage = 360°
T1 = 0;
T2 = 360;
T = linspace(T1,T2,360);
M = T-T1; 

%-------------------------------------------------------------------------
%%
% Beginn Rechnung

% Summation Besselreihe k= 1...100;
jend =6;
kend = 100;
for j =1:jend
    Summe = 0;
    ex(j) = 0.15*(j);
    for k = 1:kend
        Summe = 2.*sind(k*M)*besselj(k, k*ex(j))/k + Summe; 
    end
    E(j,:) = Summe(:) + M(:);
    upsilon(j,:)= 2*atand(sqrt((1+ex(j))/(1-ex(j)))*tand(E(j,:)/2));
    x(j,:) = (1-ex(j).*ex(j)).*cosd(upsilon(j,:))./(1+ex(j)...
             .*cosd(upsilon(j,:)));
    y(j,:) = (1-ex(j).*ex(j)).*sind(upsilon(j,:))./(1+ex(j)...
             .*cosd(upsilon(j,:)));
end

% -------------------------------------------------------------------------
% Graphische Ausgabe

header1='Exzentrische Anomalie - Besselreihe'; 
figure('Name',header1);
hold on;
for j=1:jend
    plot(T,E(j,:)-M,'Color',Colors(j,:),'LineWidth',2);
    lgdstr(j,:)=strcat('{\it e} = ',sprintf('% 4.2f \n ', ex(j)));
end
legend(lgdstr);
XTickMode = 'manual';
xticks([0 60 120 180 240 300 360])
xlim([0 360]); 
ylim([-1 1]);
grid on;
header2='Abweichung{\it E} von{\it M}';
title(header2);
xlabel('Grad')
ylabel('Abweichung in Grad');
legend boxoff;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

header1='Bahnkurven - Besselreihe'; 
figure('Name',header1);
hold on;
for j=1:jend
    plot(x(j,:),y(j,:),'Color',Colors(j,:),'LineWidth',2);
    lgdstr(j,:)=strcat('{\it e} = ',sprintf('% 4.2f \n ', ex(j)));
end
SonneFP;
legend(lgdstr);
axis equal;
xlim([-2 1]); 
ylim([-1 1]); 
grid on;
header2='Bahnkurven nach Besselreihe';
title(header2);
xlabel('x ')
ylabel('y ');
legend boxoff;
set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);

% Vergleich der Genauigkeit für verschiedene Exzentrizitäten
header1='Genauigkeit der Besselreihe (Summe bis n)'; 
figure('Name',header1);
nend = [2 4 6 10 50 1000]; % Summation 1...nend , 1000 entspricht unendlich
jend2 =6;
% lgdstr2 = strings([jend2,15]);

for iPlot = 1:2
    exz = 0.2*iPlot^2;
    subplot(1,2,iPlot);
    for j =1:jend2
        Summe = 0;
        for k = 1:nend(j)
            Summe = 2.*sind(k*M)*besselj(k, k*exz)/k + Summe; 
        end
        E2(j,:) = Summe(:);
    end
    hold on;
    for j=1:jend2
        plot(T,E2(j,:)-E2(jend2,:),'Color',Colors(j,:),'LineWidth',2);
        lgdstr2(j,:)=strjoin(["{\it n} = ", num2str(nend(j))]);
    end
    legend(lgdstr2,'location', 'southeast');
    XTickMode = 'manual';
    xticks([0 60 120 180 240 300 360])
    xlim([0 360]); 
    ylim([-0.005*100^(iPlot-1) 0.005*100^(iPlot-1)]);
    grid on;
    header2=sprintf('% 3.1f ',exz);
    header2=strcat('Exzentrizität {\it e} = ',header2);
    title(header2);
    xlabel('Grad')
    ylabel('Abweichung in Grad');
    legend boxoff;
    set(gca,'XGrid','on', 'YGrid', 'on', 'Fontsize', 18, 'linewidth', 1);
end

% Man sieht sehr schön, dass für kleine Exzentrizitäten die Besselreiche
% sehr schnell (nend >= 5) konvergiert, bei großen Exzentrizitäten muss man
% schon bis nend >= 50 aufaddieren.

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------