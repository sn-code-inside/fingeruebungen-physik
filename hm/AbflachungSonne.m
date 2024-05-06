% -------------------------------------------------------------------------
% AbflachungSonne.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Abflachung der Sonne bei Sonnenuntergang als Funktion der Hoehe.
% Courtesy of Bernd Geh.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
 
% % https://www.timeanddate.com/weather/usa/kailua-kona/historic?month=11&year=2016
% P=101.59; %Pressure in kPa (@Time of picture)
% T=27; % Temperature in degrees (@Time of picture)
% OH=25; % Obeserver Elevation in m
% R=6378000; % Earth radius in m
h=linspace(1.25,-1.75,6); %Sun elevation in degrees

% Sun outline
phi=linspace(0,2*pi,200);
Dsun = 32.3/60;             % Sun radius in deg
cphi=0.5*Dsun*cos(phi);     % x-xoordinate of sun
sphi=0.5*Dsun*sin(phi);     % y-xoordinate of sun

figure;
% One or 2 Pictures
pic_flag=0; % 0: Only one Picture 1: Two Pictures
if pic_flag
    rgb1 = imread('SunsetHawaii.jpg');
    rgb2 = imread('SunsetHawaii01.jpg');
    rgb=rgb1/2+rgb2/2; % Loading both sun pictures
else
    rgb = imread('SunsetHawaii.jpg'); % Only one sun picture
end

pixres = 0.00101; % Calibrated degrees per Pixel  
xra=pixres*size(rgb,1)/2;
yra=pixres*size(rgb,2)/2;
imagesc([-xra,+xra],[-yra,yra]+0.605,flip(rgb,1),'AlphaData',0.6); 
axis xy; 

hold on;
grid on;
sh = 0.6 ;   % Shift for plot
for k=1:length(h)
  cys=h(k)+sphi;
  R1=1.02./tand(cys+10.3./(cys+5.11))/60; %Refraction formula 3.76
  cy=cys+R1;
  if k < length(h)-1 
      plot(cphi-sh,cy,'Color',Colors(k,:),'LineWidth',3,'LineStyle',':'); 
      patch(cphi-sh,cy,Colors(10,:),'FaceAlpha',.75);  
  end
  plot(cphi+sh,cys,'Color',Colors(k,:),'LineWidth',2); 
end

% Location 20.0625N, 155.8498W, Nov 10, 2016, ~5:42:32pm 
h0= -0.3492; % Elevation angle input (according location and time)
cys0 = h0 +sphi;
R1=1.02./tand(cys0+10.3./(cys0+5.11))/60; %Refraction formula 3.76
cy0=cys0+R1;
plot(cphi,cy0,'Color','k','LineWidth',2,'LineStyle',':');  

if pic_flag==1 %Plot second sun
    % Location 20.0625N, 155.8498W, Nov 10, 2016, ~5:38:07pm earlier
    h0=  0.633; % Elevation angle input (according location and time)
    cys0 = h0 +sphi;
    R1=1.02./tand(cys0+10.3./(cys0+5.11))/60; %Refraction formula 3.76
    cy0=cys0+R1;
    plot(cphi,cy0,'Color','k','LineWidth',2,'LineStyle',':');
end

% Horizon angle (to match the picture) 
hor = -0.0907;
plot([-1,1],[hor,hor],'Color',Colors(3,:),'LineWidth',2,'LineStyle','-');

xh = [-1 1 1 -1];
yh = [-1 -1 hor hor];
patch(xh,yh,Colors(8,:),'FaceAlpha',.65);  

axis square;
ylabel('Höhe in °');
xlabel('Azimut (relativ) in °');
xlim([-0.88,0.88]);
ylim([-0.27,1.49]);
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------