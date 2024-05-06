% -------------------------------------------------------------------------
% MarsOpposition.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Marsoppositionen im Zeitraum 2000-2020. 
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Initialisierung
AE=149;
pi2=2*pi;

dt0 = datetime('2000-01-01 12:00:00');
T0 = juliandate(dt0);
Aequi = 'J2000';
[BaPa,BaPadot]=OrbitParameter(T0, Aequi);

dt1 = datetime('2000-01-01 00:00:00');
T1 = juliandate(dt1); % Julianisches Datum
MJuDa1 = juliandate(dt1,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr1 =  string(dt1,'dd.MM.yyyy');
Year1     =  string(dt1,'yyyy');
dt2 = datetime('2020-01-01 00:00:00');
T2 = juliandate(dt2); % Julianisches Datum  %Bedeckungsdatum
MJuDa2 = juliandate(dt2,'modifiedjuliandate');% Modif. Julianisches Datum
DatumStr2 =  string(dt2,'dd.MM.yyyy');
Year2     =  string(dt2,'yyyy');


%-------------------------------------------------------------------------
%Begin Rechnung
maxP=10000;
T_vector=linspace(T1,T2,maxP);

%Berechnung Erdkreis
for k=1:12  
    dt3 = datetime([2020, k, 1, 0, 0,0]);
    MJuDa3 = juliandate(dt3,'modifiedjuliandate');% Modif. Julianisches Datum
    DateStr3 =  string(dt1,'dd.MM.yyyy');
    DatumStr3(k).month = DateStr3 ;
    Monat(k)=juliandate(dt3); % Julianisches Datum
end

Erde=PlanetPQR(Monat, BaPa, BaPadot, 3);

% Berechnung Planetenposition als Funktion der Zeit nach Keplerlösung
for k=3:4
    Planets(k)=PlanetPQR(T_vector, BaPa, BaPadot, k);
end
for iPlot = 3:4
    iPn=deg2rad(BaPa.iP(iPlot));
    PiqPn=deg2rad(BaPa.PiqP(iPlot));
    OmegaPn=deg2rad(BaPa.OmegaP(iPlot));
    omegaP=PiqPn-OmegaPn;
    Planets(iPlot).per=[BaPa.aP(iPlot)*(1-BaPa.eP(iPlot));0;0];
    Planets(iPlot).per=mtimes(PQR(-OmegaPn,-iPn,-omegaP),Planets(iPlot).per);
end

Planets(3).ekl(2,:)= wrapTo2Pi(Planets(3).ekl(2,:));
Planets(4).ekl(2,:)= wrapTo2Pi(Planets(4).ekl(2,:));

%Berechnung Erdkreis
for k=1:12
    dt3 = datetime([2020, k, 1, 0, 0,0]);
    MJuDa3 = juliandate(dt3,'modifiedjuliandate');% Modif. Julianisches Datum
    DateStr3 =  string(dt3,'MMM');
    DatumStr3(k).month = DateStr3 ;
    Monat(k)=juliandate(dt3); % Julianisches Datum
end
Erde=PlanetPQR(Monat, BaPa, BaPadot, 3);

%Suche der Oppositionspositionen und Abstand
found = 0;
for k=2:maxP
      if (Planets(4).ekl(2,k) - Planets(4).ekl(2,k-1)) < 0 
        fM =pi2; 
      else
        fM=0;
      end  
      if (Planets(3).ekl(2,k) - Planets(3).ekl(2,k-1)) < 0 
        fE =pi2;
      else
        fE=0;
      end  
      if (Planets(4).ekl(2,k) + fM - Planets(3).ekl(2,k) - fE)*(Planets(4).ekl(2,k-1)-Planets(3).ekl(2,k-1)) < 0 
          found = found +1;
          opp_x(found)= Planets(4).xyz(1,k);
          opp_y(found)= Planets(4).xyz(2,k);
          DateFound = datetime(Planets(4).Time(k),'convertfrom','juliandate');
          opp_t(found,:) = string(DateFound,'dd-MM-yyyy');
          findex(found)=k;
          abstand(found)=sqrt((Planets(4).xyz(1,k)-Planets(3).xyz(1,k))^2+(Planets(4).xyz(2,k)-Planets(3).xyz(2,k))^2+(Planets(4).xyz(3,k)-Planets(3).xyz(3,k))^2);
          abstand(found)=abstand(found)*AE;
          dist(found).str =sprintf('%3d ',round(abstand(found)));
          dist(1).str =sprintf('%3d Mio km',round(abstand(1)));
          delta(found)=1;
          alt(found)  =0.5;
      end
    foundmax=found;
 end
dist(round(foundmax/2)).str =sprintf('%3d Mio km',round(abstand(round(foundmax/2))));


%------------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphische Ausgabe

%Berechnung Frühlingspunkt

header1='Bahnen von Mars und Erde auf die Ekliptik projiziert';
figure('Name',header1);

%Ausgabe Erdkreis
for k=1:12
    p=plot(Erde.xyz(1,k),Erde.xyz(2,k),'+','Color', Colors(3,:));
    text(Erde.xyz(1,k),Erde.xyz(2,k), DatumStr3(k).month ,'Color', Colors(3,:),'FontSize',16);
    hold on;
end
for iPlot = 3:4
    plot(Planets(iPlot).xyz(1,:),Planets(iPlot).xyz(2,:),'Color', Colors(iPlot,:),'LineWidth',1);
    hold on
    axis equal
end
SonneFP;
%Ausgabe Oppositionen 
for k=1:foundmax
    x1 = Planets(3).xyz(1,findex(k));
    x2 = Planets(4).xyz(1,findex(k));
    y1 = Planets(3).xyz(2,findex(k));
    y2 = Planets(4).xyz(2,findex(k));
    x=linspace(x1,x2,100);
    y=linspace(y1,y2,100);
    plot(x,y,'Color', Colors(4,:));
    text((x1+x2)/2, (y1+y2)/2,dist(k).str,'Color', Colors(4,:),'FontSize',14);
    text(Planets(4).xyz(1,findex(k)),Planets(4).xyz(2,findex(k)),opp_t(k,:),'Color', 'k','FontSize',14);
end
%Positionen Perihel
for iPlot = 3:4
    plot(Planets(iPlot).per(1),Planets(iPlot).per(2),'o','Color', Colors(iPlot,:),'LineWidth',2);
    text(Planets(iPlot).per(1)+0.15,Planets(iPlot).per(2)+0.1,'P','HorizontalAlignment','right','Color', Colors(iPlot,:),'FontSize',14);
end

header2='Marsoppositionen ' + Year1 + ' - ' + Year2;
title(header2,'FontSize',12);
ylim([-2 2])
xlim([-2 2]);
grid on;
grid minor,
xlabel('x in AE')
ylabel('y in AE');
set(gca,'FontSize',16);


%%
%-----------------------------------------------------------------
%Ausgabe Bildschirm
fprintf('\n  Datum   \t l[°] \t\t b[°]\t\t dr Mio km \t aM Mio km \t beta[°] \t delta[°] \t h[°] \t Bildwinkel["]');
eps0=23.43929111;
phi =48.8;
RM  =6772;

for k=1:foundmax
    out1 =  rad2deg(Planets(4).ekl(2,findex(k)));%helioz. Länge Mars
    out2 =  rad2deg(Planets(4).ekl(3,findex(k)));%helioz. Breite Mars
    out3 =  abstand(k);%abstand Erde Mars
    out4 =  Planets(4).ekl(1,findex(k))*AE;%Marsabstand Sonne
    out5 =  rad2deg(asin(Planets(4).xyz(3,findex(k))*AE/out3));%beta
    out6 =  asind(sind(out5)*cosd(eps0)+ cosd(out5)*sind(out1)*sind(eps0));%delta
    out7 =  asind(cosd(out6)*cosd(phi)+sind(out6)*sind(phi)); %Alt
    out8 =  atand(RM*1e-6/out3)*3600; %Bildwinkel in Bogensekunden
    fprintf('\n%s \t %03.2f  \t %+2.2f   \t %+3.1f  \t %+3.1f  \t %+2.2f   \t %+2.2f  \t %+2.2f \t %+2.2f ', opp_t(k,:), out1, out2, out3, out4, out5, out6, out7, out8);   
end
fprintf('\n  \n');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------
