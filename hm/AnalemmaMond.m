% -------------------------------------------------------------------------
% AnalemmaMond.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial wurde von Kenneth von Buenau und Michael
% Kaschke erstellt.
% -------------------------------------------------------------------------
% Analemma des Mondes von der Erde 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung Epoche, Auswahl verschiedener Zeitpunkte, die
% unterschiedliche Analemmaformen ergeben
dt2 = datetime('2020-03-26 18:00:00');
dt1 = datetime('2014-03-15 18:00:00');
% dt1 = datetime('today');
T1  = juliandate(dt1);
t1  = Jd2JJht(T1);

epsE    = EpsErde(T1);  
epsErad = deg2rad(epsE); 
jend = 2;
NrJDpOrbit = 27.32166; % Anzahl JD pro siderischen Mondmonat

% Datumsvektor fuer Analemma (ca. 2 Monate)
NrDays         = 30;   % Julianische Tage für die Berechnung
PointsperDay   = 96;   % Anzahl Berechnungspunkte
NrPoints       = NrDays*PointsperDay;
for k=1:NrDays
   MondTage(k) = 1+(k-1)*PointsperDay;
end
for j=1:jend
    if j==1 
        T1 = juliandate(dt1);
    else
        T1 = juliandate(dt2);
     end
    t1  = Jd2JJht(T1);
    TA          = T1;
    TF          = T1 + NrDays;
    Tv          = linspace(TA,TF,NrPoints);   % Jul. Tage
    tv          = Jd2JJht(Tv);              % Julian. Jahrhunderte 
    dtv(j,1:length(Tv)) = datetime(Tv,'convertfrom','juliandate');
    time = (Tv-TA)/NrJDpOrbit;

    for tag = 1:NrPoints
        MoonPos = MondExakt(tv(tag),epsErad,'rE');
        pos(tag,:) = MoonPos.equ;           % Mondposition in geoz-aequ. Koord.
        meanmoonpos(tag) = 2*pi*Frac(time(tag));
                                            % Definition eines mittl. Mondes
    end
    zgl = -(pos(:,2)-meanmoonpos(:))/pi*720; % Unterschied der Rektaszensionen
    Dec(j,:) = rad2deg(pos(:,3));                 % Deklination des Mondes

    for iz = 1:jend
        zgl = wrapTo180(zgl/720*180)/180*720;
        % Durch Integration ermittelter Mittelwert der ZGL
        integral = trapz(Tv,zgl);
        offset = integral/(TF-TA);
        ZGL(j,:) = zgl-offset;
    end
end

% Vollmondtage 2014
VollmondStr = ...
['16-Mar-2014 17:08:11'; ...    
 '08-Apr-2020 02:34:44'];    
VollMond     = datetime(VollmondStr,'InputFormat','dd-MMM-yyyy HHH:mm:ss');
VollMondTag  = day(VollMond,'dayofmonth');
VollMondMonat= month(VollMond,'monthofyear');
for index = 1:length(VollMond)  
    for tag = 1: length(dtv(index,:))
        MondTag  = day(dtv(index,tag),'dayofmonth');
        MondMonat  = month(dtv(index,tag),'monthofyear');
        if MondTag == VollMondTag(index) && ...
             MondMonat == VollMondMonat(index)
           IstVollMond(index) = tag;
        end
    end
end



% Bild Analemma
figure();
plot(ZGL(1,:),Dec(1,:),'-+','MarkerIndices', MondTage,'Color', ...
     Colors(4,:),'LineStyle','-','LineWidth',2);  
 %Simulation
hold on
plot(ZGL(2,:),Dec(2,:),':d','MarkerIndices', MondTage,'Color', ...
     Colors(10,:),'LineStyle',':','LineWidth',2);  
grid on
plot(ZGL(1,IstVollMond(1)), Dec(1,IstVollMond(1)),'o','MarkerSize',10,...
    'LineWidth',2,'Color',Colors(2,:));
plot(ZGL(2,IstVollMond(2)), Dec(2,IstVollMond(2)),'o','MarkerSize',10,...
    'LineWidth',2,'Color',Colors(2,:));
ylim([-40 40]);
xlim([-40 40]);
title('Mond von der Erde gegen mittleren Mond');
ylabel('Deklination')
xlabel('ZGL in min');
% for k=1:5:length(zgl)
% hold on
%  p(2) = plot(ZGL(1,k),Dec(1,k));
%  p(2).Color = Colors(4,:);
%  p(2).LineWidth = 3;
%  p(2).Marker = 'o';
%  p(2).Visible = 'on';
%  pause(0.002)
%  p(2).Visible = 'off';
% end
% for k=1:5:length(zgl)
% hold on
%  p(2) = plot(ZGL(2,k),Dec(2,k));
%  p(2).Color = Colors(10,:);
%  p(2).LineWidth = 3;
%  p(2).Marker = 'o';
%  p(2).Visible = 'on';
%  pause(0.002)
%  p(2).Visible = 'off';
% end
legend(strcat(string(dtv(1,1),'dd.MM.yyyy'),' - ',...
        string(dtv(1,length(dtv(1,:))),' dd.MM.yyyy')),...
       strcat(string(dtv(2,1),'dd.MM.yyyy'),' - ',...
        string(dtv(2,length(dtv(2,:))),' dd.MM.yyyy')),... 
        'Vollmondposition', 'location','southwest');
legend box off;
set(gca,'Fontsize',16);


% %%
% %Laengerer Zeitraum um Präzession der Apsidenlinie zu zeigen
% %Datumsvektor fuer Zeitbereich (zwei Jahre)
% NrDays = 10*365;   % Julianische Tage für die Berechnung
% PointsperDay = 10*365;   % Anzahl Berechnungspunkte
% TF           = T1 + NrDays;
% Tv2          = linspace(T1,TF,PointsperDay);   %Jul. Tage
% tv2          = Jd2JJht(Tv2);              %Julian. Jahrhunderte 
% dtv2 = datetime(Tv2,'convertfrom','juliandate');
% time = (Tv-T1)/NrJDpOrbit;
% 
% for tag = 1:PointsperDay
%     MoonPos2 = MondExakt(tv2(tag),epsErad,'rE');
%     pos2(tag,:) = MoonPos2.equ;           %Mondposition in geoz-aequ. Koord.
% end
% Dec2 = rad2deg(pos2(:,3));                %Deklination des Mondes
% 
% % Bild Deklination
% figure();
% plot(dtv2, Dec2,'Color', Colors(5,:),'LineStyle','-','LineWidth',1);  
% grid on,
% ymin =-40;
% ymax =40;
% ylim([ymin ymax]);
% title('Deklination des Monds von der Erde');
% legend(strcat(string(dtv2(1),'dd.MM.yyyy'),' - ',...
%        string(dtv2(length(dtv2)),' dd.MM.yyyy')),'location','southwest');
% legend box off;
% ylabel('Deklination')
% xlabel('Datum');
% set(gca,'Fontsize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
