% -------------------------------------------------------------------------
% Sonnenfinsternis01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Das vorliegende Programmaterial basiert in Teilen 
% auf C++ Programmstrukturen und der Beschreibung/Erlaeuterung in  
%
% "Astronomie mit dem Personal Computer"
%
% von Oliver Montenbruck und Thomas Pfleger (Springer 1999). 
% Genehmigung des Verlages und der Autoren liegt vor.
% -------------------------------------------------------------------------
% Berechnet die Zeiten des Neu/Vollmondes in einem bestimmten Zeitraum und
% bestimmt, ob zu der Neumond/Vollmondzeit eine SoFi/MoFi stattfinden kann.
% Fuer die Koordinaten Sonne und Mond werden genauere Reihen 
% verwendet.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
SunData = PerturbImport('SonnePos.csv');

%%
%Beginn Programm

MPhase = 0; 
% Parametrisierte Funktion zur Berechnung der relativen Phasenlage der 
% Laengen von Mond und Sonne
myfunc = @(x,a,b) PhaseFunc(x,SunData,MPhase);  
a  = SunData; b  = MPhase; % Parameter
% Function Handle
func = @(T) myfunc(T,a,b);    % Funktion von T mit a, b Parameter

Jahr0 = 2021;
n_years  = 30;  % Anzahl der Jahre
n_months = 13;  % Anzahl der Monate 
delt = 7/36525; % Wochenweise Berechnung
tzero = zeros();
 
titlestr = strings(2,25);
titlestr(1,:) = 'Sonnenfinsternisse';
titlestr(2,:) = 'Mondfinsternisse';
% Bestimmung des Auftretens von Finsternissen in verschiedenen Jahren

for iPlot = 1:2
    disp('berechnet die Daten');
    MPhase= (iPlot-1)*2; %0 Neumond  2 Vollmond
   % Zunaechst Abfrage der Zeitpunkte der Neumond/Vollmondzeitpunkte
    for k= 1:n_years
        Jahr = Jahr0+k;
        dt0 = datetime(Jahr,12,01); %Startdatum ET
        T0  =  juliandate(dt0);
        eps0 = deg2rad(EpsErde(T0));
        t0 = Jd2JJht(T0); %Zeit in Jul. Jahrhunderten seit JD2000
        t1 = t0+delt;
        for j = 1:n_months
              D0 = PhaseFunc(t0,SunData,MPhase);%Berechnungszeitpkt 1.Woche  
              D1 = PhaseFunc(t1,SunData,MPhase);%Berechnungszeitpkt 2.Woche 
              while (D0*D1 > 0) || (D1<D0) %kein Neumond
                  t0 = t1;
                  D0 = D1;
                  t1 = t1+delt;
                  D1 = PhaseFunc(t1,SunData,MPhase);
              end
              tber =[t0,t1]; % Wochenbereich fuer Neumond/Vollmond
              Eclipse(iPlot).t(k,j) = fzero(func,tber); %genaue Zeit
              t0 = t1;
              t1 = t0 + delt;
        end
    end
    fprintf(' \n');
end

for iPlot = 1:2
    % Nun ueberpruefung auf moegliche Bedeckung der Sonne durch Neumond bzw. 
    % Mond durch Erdschatten
    % ueberpruefung der ekliptikalen Breite in EclipseFlag
    disp('berechnet die Daten'); %Warteinfo
  for k= 1:n_years
      fprintf('\n %s  \n', titlestr(iPlot));
      for j = 1:n_months
        tzero(k,j) = Eclipse(iPlot).t(k,j); 
        T = JJht2Jd(tzero(k,j));
        DelTSec = ETminusUT(T);  %Unterschied UT zu ET
        T = T - DelTSec/86400;
        Tag =   datetime(T,'ConvertFrom','juliandate');
        epsE = deg2rad(EpsErde(T));
        MoonPos = MondExakt(tzero(k,j) ,epsE,'RE'); 
        Eclipse(iPlot).Tag(k,j) = Tag;
        Eclipse(iPlot).Flag(k,j) = EclipseFlag(MoonPos.ekl(3),iPlot);
        %Ausdruck
        outstr(j,:) = string(Tag)+" UT"+' '+  ... 
                      string(Eclipse(iPlot).Flag(k,j));
        fprintf(' %s \n',outstr(j,:));
     end
 end
  fprintf('\n \n');  
end


%%
% ------------------------------------------------------------------------
%
% Graphische Ausgabe

for iPlot = 1:2
    figure();
    title(titlestr(iPlot));
    hold on;
    first_Total = true; 
    first_Partiell = true; 
    for k= 1:n_years
        for j = 1:n_months
            Tag = Eclipse(iPlot).Tag(k,j);
            if strfind(Eclipse(iPlot).Flag(k,j),"Total")               
                 if first_Total == true
                     first_Total = false;  
                     p2 = plot(day(Tag,'dayofyear'),year(Tag),...
                        'd','Color',Colors(4-2*(iPlot-1),:), ...
                        'LineWidth',2);
                 else
                      plot(day(Tag,'dayofyear'),year(Tag),...
                     'd','Color',Colors(4-2*(iPlot-1),:),'LineWidth',2);
                 end      
            else
               if strcmp(Eclipse(iPlot).Flag(k,j),"Partiell")
                  if first_Partiell == true
                      first_Partiell = false;     
                    p1 = plot(day(Tag,'dayofyear'),year(Tag),...
                        's','Color',Colors(3+2*(iPlot-1),:), ...
                        'LineWidth',2);  
                  else 
                      plot(day(Tag,'dayofyear'),year(Tag),...
                        's','Color',Colors(3+2*(iPlot-1),:),'LineWidth',2);
                  end     
               end
            end
        end
    end
    for k=1:13
        XDataTick(k) = datetime(2000,k,1);
    end
     XDataTick(13) = datetime(2001,1,1);
    ax = gca;
    datetick('x','mmm','keepticks')
    legend([p1, p2], ["Partiell", "Total"], 'location', 'southeast');
    xlim ([1 395]);
    ylim ([Jahr0 Jahr0+n_years+2]);
    xlabel('Monat');
    ylabel('Jahr');
    grid on;
    set (gca,'FontSize',16);
end
 
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------