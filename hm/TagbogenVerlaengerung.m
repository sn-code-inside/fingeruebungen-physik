% -------------------------------------------------------------------------
% TagbogenVerlaengerung.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Tagbogenverlaengerung durch Refraktion und scheinbaren 
% Sonnendurchmesser bei Sonnenaufgang. Gilt analog fuer Sonnenuntergang
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
t1 = datetime('2000-01-1 0:0:00');
% Julianisches und Modifiziertes Julianisches Datum
T = juliandate(t1);
MJuDa =juliandate(t1,'modifiedjuliandate');

%Refraktion und scheinbarer Durchmesser in °
hr=50/60;
%Umrechnung Bogenminuten in Zeitsekunden
hrz=hr*3600/15;
% Schleifen ueber 6 geografische Breiten und 365 Tage 
deltatH1 = zeros(7,366);
titlestr   =  strings(7,1);
day=1:366;
for n=1:7
    phi= 0.001+(n-1)*10;
    for m=1:366
        %Naeherungsweise Berechnung der Deklination nach Sinusformel
        delta = rad2deg(0.4095)*sin(0.016906*(m-80.086));
        tauH=acos(-tand(phi)*tand(delta));     
        if ~isreal(tauH)
            tauH=NaN;
        end   
        %Tagbogenverlaegerung
        deltatH1(n,m)=hrz/(cosd(phi)*cosd(delta)*sin(tauH));
    end
    titlestr(n)=['Breite ', sprintf('%2.0f°     ',phi)];
end

%Ausgabe

figure('Name','Berechnetes Delay zw. geometrisch-astronom. und beobacht. SA');
for k=1:6
    plot(day,deltatH1(8-k,:)/60,'Color',Colors(8-k,:),'LineWidth',1);
    hold on;
end
header='Tagbogenverlängerung durch Refraktion und scheinbaren Radius';
text(10,11,header,'FontSize',18);

xlim([1 366]);
ylim([2,12]);
xlabel('Tag');
ylabel('Minuten');
lgd=legend(titlestr(6),titlestr(5),titlestr(4),titlestr(3),titlestr(2),titlestr(1),'location','south','numcolumns',4);
lgd.FontSize=18;
legend boxoff;
grid on;
set(gca,'FontSize',18);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------