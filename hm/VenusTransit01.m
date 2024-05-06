% -------------------------------------------------------------------------
% VenusTransit01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die Koordinaten der Venus und Merkur bei Transits
% Die Zeiten müssen naeherungsweise z.B. aus InnerePlanetenPQR.m bekannt 
% sein. Die Daten werden mit JPL-Daten verglichen. 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Daten von JPL Horizons
fid=fopen('Venus20040608.dat','r');
if fid ==1 
     disp('File open not successful');
else
    VenusJPL = ImportfileJPL('./IncludeFolder/Venus20040608.dat', 2, inf);    
end
closeresult =fclose(fid);
if closeresult ==0
%     disp('Data successfully loaded');
else
    disp('Dataload not successful');
end

fid=fopen('Sonne20040608.dat','r');
if fid ==1 
    disp('File open not successful');
else
    SunJPL = ImportfileJPL('Sonne20040608.dat', 2, inf);    
end
closeresult =fclose(fid);
if closeresult ==0
%     disp('Data successfully loaded');
else
    disp('Dataload not successful');
end

%-------------------------------------------------------------------------
%Berechnung der Koordinaten z.B. am 8.6.2004 (Termin des Transits)
%Zeiten in UT
dt1 = datetime('2004-06-08 04:00:00'); %Venustransit 2004 
dt2 = datetime('2004-06-08 14:00:00');

% dt1 = datetime('2012-06-05 20:00:00'); %Venustransit 2012
% dt2 = datetime('2012-06-06 07:00:00');
% 
% dt1 = datetime('2019-11-11 12:00:00'); %Merkurtransit 2019
% dt2 = datetime('2019-11-11 19:00:00');

T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
T2 = juliandate(dt2); % Julianisches Datum  %Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');

Aequi = 'Datum';
Aequi = 'J2000';
% Einlesen der Bahnparameter und deren Ableitung 
[BaPa,BaPadot]=OrbitParameter(T2, Aequi);

if year(dt1) == 2004 
    compJPL = 1;
else
    compJPL = 0;
end

%-------------------------------------------------------------------------
%Bestimmung der Schrittweite zur Zeitberechnung
St_p_H=4;
if day(dt1) == day(dt2)
    Nhour=hour(dt2)-hour(dt1);
else
    Nhour = 24-hour(dt1) + hour(dt2);
end
maxP=Nhour*St_p_H+1; 
Lab_p_H=1;  % jede Stunde ein Tick setzen auf der Bahn
NrLabSteps = St_p_H/Lab_p_H;
NrLabels=Nhour*Lab_p_H+1;

T_vector=linspace(T1,T2,maxP);
MJD_vector=linspace(MJuDa1,MJuDa2,maxP);
Stunden=linspace(hour(dt1),hour(dt1)+Nhour,NrLabels);
DelTSec = ETminusUT(T1)/86400;

%__________________________________________________________________________
% Berechnung der Koordinaten nach PQR nach Keplerloesung
for k=1:3
    temp(k)=PlanetPQR(T_vector+DelTSec, BaPa, BaPadot, k);
end

%Optional: Ausgabe der ekliptikalen Koordinaten im Command Window
for k=1:maxP
%     DateObs = datetime(MJD_vector(k),'ConvertFrom','modifiedjuliandate');
%     Datestr = string(DateObs, 'dd-MM-yyyy HH:mm ');
%     l1=wrapTo360(rad2deg(Planets(1).ekl(2,k)));
%     b1=rad2deg(Planets(1).ekl(3,k));
%     l2=wrapTo360(rad2deg(Planets(2).ekl(2,k)));
%     b2=rad2deg(Planets(2).ekl(3,k));
%     l3=wrapTo360(rad2deg(Planets(3).ekl(2,k)));
%     b3=rad2deg(Planets(3).ekl(3,k));
%     fspec="\n %S lM %8.6f lV %8.6flE %8.6f bM %8.6f bV %8.6f bE %8.6f";
%     fprintf(fspec, Datestr, l1,l2,l3,b1,b2,b3);
end

% Umrechnung in aequatoriale Koordinaten 
epsErad = deg2rad(EpsErde(T2));
Erde = temp(3);
for k=1:3
    Planets(k) = Convert2Equ(temp(k), Erde, epsErad);
end

%__________________________________________________________________________
% Berechnung und Ausgabe der relativen Koordinaten vor der Sonnenscheibe
% in Rektaszension und Deklination
for k=1:maxP
    DateObs = datetime(MJD_vector(k),'ConvertFrom','modifiedjuliandate');
    Datestr = string(DateObs, 'dd-MM-yyyy HH:mm ');
    %Umrechnung in relative aequatoriale Koordinaten (Aequinoktium J2000)
    RAV(k) = rad2deg(Planets(2).equ(2,k) - Planets(3).equ(2,k));
    DekV(k) = rad2deg(Planets(2).equ(3,k) - Planets(3).equ(3,k));
    RAM(k) = rad2deg(Planets(1).equ(2,k) - Planets(3).equ(2,k));
    DekM(k) = rad2deg(Planets(1).equ(3,k) - Planets(3).equ(3,k));
      
    if compJPL  == 1
        RAV_JPL(k) = VenusJPL.Alpha(k) - SunJPL.Alpha(k);
        DekV_JPL(k) = VenusJPL.Delta(k) - SunJPL.Delta(k);
    end
%     aM = StrHMS(-RAM(k)/15);
%     dM = StrDMS(DekM(k));
%     aV = StrHMS(-RAV(k)/15);
%     dV = StrDMS(DekV(k));
%     fprintf('\n %s  aMerkur: %s  dMerkur%s  aVenus: %s  dVenus %s ', ...
%             Datestr, aM, dM ,aV, dV);
end
figure();
plot(-RAV, DekV,'LineWidth',2,'Color',Colors(2,:));
hold on;
if compJPL == 1 
    plot(-RAV_JPL, DekV_JPL,'LineWidth',2,'Color',Colors(3,:));
end
plot(-RAM, DekM,'LineWidth',2,'Color',Colors(1,:));
PlotCircle(0,0,0.27, Colors(10,:),3);
plot(-RAM, DekM,'-+','MarkerIndices',1:NrLabSteps:length(DekM),...
      'LineWidth',1,'Color',Colors(1,:));
plot(-RAV, DekV,'-+','MarkerIndices',1:NrLabSteps:length(DekV),...
      'LineWidth',1,'Color',Colors(2,:));
if compJPL == 1 
    plot(-RAV_JPL, DekV_JPL,'-+','MarkerIndices',...
         1:NrLabSteps:length(DekV),'LineWidth',1,'Color',Colors(3,:));
end
if compJPL == 1
    lgd=legend('Venus Bahn','JPL','Merkur Bahn','location','northeast');
else
    lgd=legend('Venus Bahn','Merkur Bahn','location','northeast');
end
legend boxoff;
grid on;
grid minor,
for k=1:NrLabels
    mylabels(k,:)=string(datetime(T_vector((k-1)*NrLabSteps +1)+eps,...
                  'convertfrom','juliandate'),'HH');
    RAVL(k)=RAV((k-1)*NrLabSteps +1);
    DekVL(k)=DekV((k-1)*NrLabSteps +1); 
    RAML(k)=RAM((k-1)*NrLabSteps +1);
    DekML(k)=DekM((k-1)*NrLabSteps +1);
 end
    h=LabelPoints(-RAVL, DekVL,mylabels,'S',0.02,1,...
           'FontSize',12,'Color',Colors(2,:));
    h=LabelPoints(-RAML, DekML,mylabels,'S',0.02,1,...
           'FontSize',12,'Color',Colors(1,:));
if compJPL ==1 
    for k=1:NrLabels
        mylabels(k,:)=sprintf('%02u', Stunden(k));
        RAVL_JPL(k)=RAV_JPL((k-1)*NrLabSteps +1);
        DekVL_JPL(k)=DekV_JPL((k-1)*NrLabSteps +1);  
    end
    h=LabelPoints(-RAVL_JPL, DekVL_JPL,mylabels,'N',0.05,1,...
                  'FontSize',12,'Color',Colors(3,:));
end
header2 = strjoin(['Venus/Merkur Transit am ',DatumStr1," (UT)"]);
title(header2,'FontSize',12);
lgd=xlabel('{\alpha} in °');
lgd.FontSize=14;
lgd=ylabel('{\delta} in °');
lgd.FontSize=14;
ylim([-0.5 0.5]);
xlim([-0.5 0.5]);
axis square;
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
