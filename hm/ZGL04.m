% -------------------------------------------------------------------------
% ZGL04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Es wird die ZGL und das Analemma für die Epochen J1246, J2000 und
% J8000 berechnet.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;

% Import data from text file
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = ";";
opts.VariableNames = ["VarName1", "eps", "ex", "varpi"];
opts.VariableTypes = ["double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
ParaLT = readtable("ParaLangzeit.csv", opts)
clear opts

% Initialisierung
pi2 = 2*pi;
dt(1) = datetime('1246-01-01 12:00:00');
dt(2) = datetime('2000-01-01 12:00:00');
dt(3) = datetime('8000-01-01 12:00:00');
DeltaT = [640,64,120000];
DeltaT = DeltaT/86400;
j=1:3;
TUT1 = juliandate(dt);
T1 = TUT1-1+DeltaT;
t1 = (T1-2451545.0)/36525; %Julianische Jahrhunderte
t1(2)=0;
MJuDa = juliandate(dt,'modifiedjuliandate');
MJuDa = MJuDa+DeltaT/86400;
titlestr   =  strings([3,25]);
eps = zeros(3,1);
ex = zeros(3,1);
varpi = zeros(3,1);
for j=1:3
    titlestr(j,:) = datestr(dt(j),'yyyy');
end
for j=1:3
    eps(j)=ParaLT.eps(1);
    for k=2:7
        eps(j) = ParaLT.eps(k)*t1(j)^(k-1)/3600+eps(j);
    end
    varpi(j)=ParaLT.varpi(1)+180;
    for k=2:7
        varpi(j) = ParaLT.varpi(k)*t1(j)^(k-1)+varpi(j);
    end
    ex(j)=ParaLT.ex(1);
    for k=2:7
        ex(j) = ParaLT.ex(k)*t1(j)^(k-1)+ex(j);
    end
    alpha(j,1)=-2*ex(j)*cosd(varpi(j))-2*ex(j)*cosd(varpi(j))*(tand(eps(j)/2))^2;
    alpha(j,2)=-5/4*ex(j)*ex(j)*cosd(2*varpi(j))+(tand(eps(j)/2))^2;
    alpha(j,3)=2*ex(j)*cosd(varpi(j))*(tand(eps(j)/2))^2;
    alpha(j,4)=-1/2*(tand(eps(j)/2))^4;
    beta(j,1)= 2*ex(j)*sind(varpi(j))- 2*ex(j)*sind(varpi(j))*(tand(eps(j)/2))^2;
    beta(j,2)=5/4*ex(j)*ex(j)*sind(2*varpi(j));
    beta(j,3)=-2*ex(j)*sind(varpi(j))*(tand(eps(j)/2))^2;
    beta(j,4)=0;
    %Für große epochale Abstände muss man die Entwicklung der mittleren
    %Länge etwas genauer approximieren. Die Formel ist bezogen auf das
    %Äquinoktium des Datums.
    L0S(j) = 180+100.46664567+36000.76982779*t1(j)+0.0003032028*t1(j)*t1(j)+t1(j)*t1(j)*t1(j)/49931000-t1(j)*t1(j)*t1(j)*t1(j)/153000000;
end

%-------------------------------------------------------------------------
%Beginn Rechnung

% ZGL1 nach Formel (3.83) 
% Achtung, da wir die veränderten Bahnparameter schon oben zeitabhängig berechnet haben,
% hier nur noch innerhalb eines Jahres (mit den entsprechenden Parametern rechnen)
Timey = 1:366;
ZGL1=zeros(3,366);
for j=1:3
  Temp = 0;
  for k=1:4
    gamma     =-24*60/pi2*sqrt(alpha(j,k)*alpha(j,k)+beta(j,k)*beta(j,k));
    f         = k*35999.372/36525; 
    g         = mod(atan2d(beta(j,k),alpha(j,k))+k*L0S(j),360)-180;
    Temp      = gamma*sind(f*Timey + g) +Temp;
  end
  ZGL1(j,:)=Temp(:);
end

% ZGL nach Formel Smart1986
% ist einfacher, da ZGL einfach als die Länge der mittleren Sonne 
% minus die Rektaszension der wahren Sonne ist
ZGL2=zeros(3,366);
for j=1:3 
    tA = (T1(j)-2451545.0)/36525;
    tE = (T1(j)+366-2451545.0)/36525;
    TV = linspace(tA,tE,366);
    L0 = pi2*Frac(280.46664567/360 +(36000.76982779*TV+0.0003032028*TV.*TV+TV.*TV.*TV/4993000-TV.*TV.*TV.*TV/153000000-TV.*TV.*TV.*TV.*TV/2e11)/360);
    M  = pi2*Frac((357.52772+35999.050340*TV+0.0001603*TV.*TV-TV.*TV.*TV/300000)/360);
    y  = (tand(eps(j)/2))^2;
    ZGL2(j,:)= y*sin(2*L0)-2*ex(j).*sin(M)+4*y*ex(j).*sin(M).*cos(2*L0)-0.5*y.*y.*sin(4*L0)-5*ex(j).*ex(j).*sin(2*M)/4; 
    ZGL2(j,:)= 24*60*ZGL2(j,:)/pi2;
end

% ZGL/Analemma numerisch Keplerlösung
for j=1:3
    MJuDa_vector = MJuDa(j)+1:MJuDa(j)+366;
    tA = T1(j)+1;
    tE = T1(j)+366;
    TV = linspace(tA,tE,366);
    % ZGL über Kepler-Lösung (Mittelpunktsgleichung) Variation eps
    % Berechnung nach Keplerlösung
    [L0, RAV,DecV]=KeplerSonneLang(TV, deg2rad(eps(j)), ex(j));
    % Berechnung des Stundenwinkels
    tau0 = 4*wrapTo180(rad2deg(L0 - RAV));
    ZGL3(j,:)=tau0;
    Deklination(j,:)= rad2deg(DecV);
    RA(j,:)=rad2deg(RAV);
end


%------------------------------------------------------------------------------
% Graphische Ausgabe
Nulllinie=zeros(366);
BeginnMonat=[ 1 32 60 91 121 152 182 213 244 274 305 335];

%Analemma
header1='Analemma'; 
figure('Name',header1);
hold on;
plot(ZGL3(1,:),Deklination(1,:),'-+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(2,:),'LineStyle','-');
plot(ZGL3(2,:),Deklination(2,:),'-+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(3,:),'LineStyle','-');
plot(ZGL3(3,:),Deklination(3,:),'-+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(4,:),'LineStyle','-');
%Näherungsformeln
% plot(ZGL1(1,:),Deklination(2,:),':+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(2,:),'LineStyle',':');
% plot(ZGL1(2,:),Deklination(2,:),':+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(3,:),'LineStyle',':');
% plot(ZGL1(3,:),Deklination(2,:),':+','MarkerIndices',BeginnMonat,'LineWidth',2,'Color',Colors(4,:),'LineStyle',':');
for j=1:3
    for k=1:12
        mylabels(k,:)=sprintf('1.%02u',k);
        xL(k)=ZGL3(j,BeginnMonat(k));
        yL(k)=Deklination(j,BeginnMonat(k));
    end
    h=LabelPoints(xL, yL ,mylabels,'S',0.05,1,'FontSize',8,'Color',Colors(1+j,:));
end
ylim([-40 30]);
xlim([-20 30]);
grid on;
xlabel('Zeit in min')
ylabel('Deklination °');
legend(titlestr(1),titlestr(2),titlestr(3),'location','southeast');
legend boxoff;
set(gca,'FontSize',18);

%Zeitgleichung
header1='Zeitgleichung';
figure('Name',header1);
% subplot(3,1,1);
plot(Timey, ZGL3(1,:), 'Color', Colors(2,:),'LineStyle','-','LineWidth',2);
hold on
plot(Timey, ZGL3(2,:), 'Color', Colors(3,:),'LineStyle','-','LineWidth',2);
plot(Timey, ZGL3(3,:), 'Color', Colors(4,:),'LineStyle','-.','LineWidth',2);
%Näherungsformeln
% plot(Timey, ZGL1(1,:), 'Color', Colors(2,:),'LineStyle',':','LineWidth',2);
% plot(Timey, ZGL1(2,:), 'Color', Colors(3,:),'LineStyle',':','LineWidth',2);
% plot(Timey, ZGL1(3,:), 'Color', Colors(4,:),'LineStyle',':','LineWidth',2);
plot(Timey, Nulllinie, 'Color', 'k','LineWidth',1);
ylim([-20 20]);
xlim([0 370]);
grid on;
xlabel('Tag')
ylabel('ZGL (Anteile) in Minuten');
lgd=legend(titlestr(1),titlestr(2),titlestr(3),'Location','south');
% lgd.FontSize =18;
legend boxoff;
set(gca,'Fontsize',18);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%%
%-------------------------------------------------------------------------
% KeplerSonneLang
% Berechnet Rektaszension und Deklination der Sonne unter Verwendung
% einer analytischen Reihenentwicklung der Keplerloesung mit relativ
% geringer Genauigkeit. Im Unterschied zu KeplerSonne.m werden hier
% verschiedene Exzentrizitaeten beruecksichtigt und es wird ein langer
% Beobachtungszeitraum gewaehlt. 
%-------------------------------------------------------------------------
%
function [L0, RAV,DecV]=KeplerSonneLang(T,eps,exz) 
% Eingabe:
%   TV       Zeit in Julianischen Datum
%   eps      Ekliptikschiefe
%   exz      Verhaeltnis der Exzentrizaet zu e = 0.016709
% Ausgabe:
%   L0       Laenge der mittleren Sonne in [rad]
%   RAV      Rektaszension der Sonne in [rad]
%   DecV     Deklination der Sonne in [rad]
%--------------------------------------------------------------------------
  pi2=2*pi;
  t=Jd2JJht(T);
% Mittlere Anomalie 
  M   = pi2*Frac((357.52772+35999.050340*t+0.0001603*t.*t-...
        t.*t.*t/300000)/360); 
% Laenge der mittleren Sonne
  L0  = pi2*Frac(280.46664567/360 +(36000.76982779*t+...
        0.0003032028*t.*t+t.*t.*t/4993000-t.*t.*t.*t/153000000-...
        t.*t.*t.*t.*t/2e11)/360);
% Mittelpunktsgleichung
  C0  = (2*exz-(exz^3)/4)*sin(M)+(5*exz*exz/4-11/24*(exz^4))*sin(2*M)+...
        13*exz^3*sin(3*M)/12;
% Ekliptikale Laenge der wahren Sonne
  L = L0+C0;
  ekl_Sun = CalcXYZfromAngles([ones(1,length(t));L;zeros(1,length(t))]);
  ekl_Sun = mtimes(R_x(-eps),ekl_Sun); % Koordinatentrafo in Equ-System
  ekl_Sun = CalcAnglesfromXYZ(ekl_Sun);
  RAV  = ekl_Sun(2,:);
  DecV = ekl_Sun(3,:); 
end
%%-------------------------------------------------------------------------
% Ende Funktion
%--------------------------------------------------------------------------

