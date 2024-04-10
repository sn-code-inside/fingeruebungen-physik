% -------------------------------------------------------------------------
% LambertProblem01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Astrodynamik" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Programm berechnet die Parameter der Ellipsenbahn 
% für das Lambert-Problem durch eine Iteration (Bisektion) am Beispiel 
% eines Flugs im ernahen Orbit.
%
% Benutzt LambertSolver1.m.

% -------------------------------------------------------------------------

% Initialisierung

clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter 
muE  = 398600.4415;   % G*ME in km^3/s^2

%% Randwertproblem Nr A:

r1     = 6800;   % in km
r2     = 6400;   % in km
gammaA = deg2rad(75) ;
TOFA   = 3000;   % in s

% Iteration des Parameters semi-latus rectum
p(1) = 2500;
p(2) = 3000;

[tFA,upsA1,upsA2,pA,exzA,ups0A]=LambertSolver1(p,TOFA,r1,r2,gammaA,muE);
V1A = sqrt(muE/pA)*[exzA*sin(upsA1), (1+exzA*cos(upsA1))];
V2A = sqrt(muE/pA)*[exzA*sin(upsA2), (1+exzA*cos(upsA2))];
v1A = vecnorm(V1A);
v2A = vecnorm(V2A);


%% Randwertproblem Nr B:

r1     = 6800;   % in km
r2     = 6400;   % in km
gammaB = deg2rad(285) ;
TOFB   = 6000;   % in s

% Iteration des Parameters semi-latus rectum
p(1) = 8000;
p(2) = 7000;
[tFB,upsB1,upsB2,pB,exzB,ups0B]=LambertSolver1(p,TOFB,r1,r2,gammaB,muE);
V1B = sqrt(muE/pB)*[exzB*sin(upsB1), (1+exzB*cos(upsB1))];
V2B = sqrt(muE/pB)*[exzB*sin(upsB2), (1+exzB*cos(upsB2))];
v1B = vecnorm(V1B);
v2B = vecnorm(V2B);
% Nachbereitung

ups = linspace(0,2*pi,361);
t   = linspace(0,2*pi,361);
r1  = r1/1000; r2  = r2/1000;  %Umrechnung in 1000 km

% Preparation
rA = pA./(1+exzA*cos(ups))/1000;
aA = pA/(1-exzA^2)/1000;
bA = aA*sqrt(1-exzA^2);
rB = pB./(1+exzB*cos(ups-ups0B))/1000;
aB = pB/(1-exzB^2)/1000;
bB = aB*sqrt(1-exzB^2);


%% Print Ausgabe

fprintf('\n')
fprintf('Randwertproblem Nr. 1:')
fprintf('\n')
fprintf('|  r1 (km)  |  r2 (km)  |  TOF (s)  | gamma (°) |  ups0 (°) |\n');
fprintf('|  %7.2f  |  %7.2f  |  %6.1f   | %7.2f   |  %7.4f  |\n',...
                r1*1000, r2*1000, TOFA, rad2deg(gammaA), rad2deg(ups0A));
fprintf('\n')
fprintf('|  a  (km)  |  p  (km)  |  ups1 (°) |  e        |  v1 (km/s)         |\n');
fprintf('|  %7.2f  |  %7.2f  |  %6.1f   |  %7.5f  | [%+7.4f,%+7.4f]  |\n',...
                aA*1000,pA, rad2deg(upsA1), exzA, V1A);
fprintf('\n')
fprintf('| v1 (km/s) | v2 (km/s) |\n');
fprintf('| %7.3f   |  %7.3f  |\n',v1A,v2A);
fprintf('\n')
fprintf('\n')
fprintf('Randwertproblem Nr. 2:')
fprintf('\n')
fprintf('|  r1 (km)  |  r2 (km)  |  TOF (s)  | gamma (°) |  ups0 (°) |\n');
fprintf('|  %7.2f  |  %7.2f  |  %6.1f   | %7.2f   |  %7.4f  |\n',...
                r1*1000, r2*1000, TOFB, rad2deg(gammaB), rad2deg(ups0B));
fprintf('\n')
fprintf('|  a  (km)  |  p  (km)  |  ups1 (°) |  e        |  v1 (km/s)         |\n');
fprintf('|  %7.2f  |  %7.2f  |  %6.1f   |  %7.5f  | [%+7.4f,%+7.4f]  |\n',...
                aB*1000,pB, rad2deg(upsB1), exzB, V1B);
fprintf('\n')
fprintf('| v1 (km/s) | v2 (km/s) |\n');
fprintf('| %7.3f   |  %7.3f  |\n',v1B,v2B);
fprintf('\n')

%% Graphik Flugellipse - Darstellung in der (r1,r2-Ebene)
[xA,yA]   = xell(aA, bA, t);
x1A = r1*cos(upsA1)+aA*exzA;
y1A = r1*sin(upsA1);
x2A = r2*cos(upsA2)+aA*exzA;
y2A = r2*sin(upsA2);
xAs = xA; yAs = yA;
xAs = flugbahn(r1, r2, exzA, aA, xAs, yAs);
% Darstellung
figure(1)
subplot(1,2,1)
plot(xA,yA,'linewidth',1,'color',Colors(2,:),'LineStyle',Style(4));
hold on
plot(xAs,yAs,'linewidth',2,'color',Colors(2,:));
plot(x1A,y1A,'s','color',Colors(2,:),...
            'MarkerSize',10,'Linewidth', 2);
plot(x2A,y2A,'d','color',Colors(2,:),...
        'MarkerSize',10,'Linewidth', 2);
% Brennpunkte
FocA1=[+aA*exzA*cos(ups0A) aA*exzA*sin(ups0A)];
FocA2=[-aA*exzA*cos(ups0A) -aA*exzA*sin(ups0A)];
PlotCircle(FocA1(1),FocA1(2),0.2,Colors(3,:),2);
PlotCircle(FocA2(1),FocA2(2),0.2,Colors(3,:),2);

%Apsidenlinie
h(1)= line([+aA*cos(ups0A) -aA*cos(ups0A)],...
     [+aA*sin(ups0A) -aA*sin(ups0A)],...
     'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(1));
%Start- und Endwert (Randwerte)
h(2)=line([FocA1(1) x1A],[FocA1(2) y1A],'color',Colors(2,:),...
          'Linewidth',2,'LineStyle',Style(2));
h(3)=line([FocA1(1) x2A],[FocA1(2) y2A],'color',Colors(2,:),...
          'Linewidth',2,'LineStyle',Style(3));
lgd =legend(h,'Apsidenlinie','r_1, t_1','r_2, t_2','location','south');
legend box off
set(lgd,'FontSize',12,'FontWeight','normal');
grid on
ylabel('y in 1000 km')
xlabel('x in 1000 km')
axis equal
xlim([-bB bB]*1.2); ylim([-aB aB]*1.2);
ttl=title('Lambert Problem');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',16);

%% Graphik Flugellipse Randwerproblem 2
[xB,yB] = xell(aB, bB, t);
% Rotation
for k=1:length(xB)
     [xB(k),yB(k)]  = rotatiere(xB(k),yB(k),ups0B);
end
x1B = r1*cos(upsB1)+aB*exzB;
y1B = r1*sin(upsB1);
x2B = r2*cos(upsB2)+aB*exzB;
y2B = r2*sin(upsB2);
[x1B,y1B]= rotatiere(x1B,y1B,ups0B);
[x2B,y2B]= rotatiere(x2B,y2B,ups0B);

xBs = xB; yBs = yB;
xBs = flugbahn(r1, r2, exzB, aB, xBs, yBs);
for k=1:length(xBs)
  if isnan(xBs(k))
      xB(k) = xB(k);
  else
      xB(k) = NaN;
  end
end
for k=1:length(xBs)
    [xBs(k),yBs(k)]= rotatiere(xBs(k),yBs(k),ups0B);
    [xB(k),yB(k)]  = rotatiere(xB(k), yB(k),ups0B);
end
subplot(1,2,2)
plot(xB,yB,'linewidth',1,'color',Colors(5,:),'LineStyle',Style(4));
hold on
plot(xBs,yBs,'linewidth',2,'color',Colors(5,:));
plot(x1B,y1B,'s','color',Colors(5,:),...
        'MarkerSize',10,'Linewidth', 2);
plot(x2B,y2B,'d','color',Colors(5,:),...
        'MarkerSize',10,'Linewidth', 2);
% Fokuspunkte 
FocB1=[+aB*exzB*cos(ups0B) aB*exzB*sin(ups0B)];
FocB2=[-aB*exzB*cos(ups0B) -aB*exzB*sin(ups0B)];
PlotCircle(FocB1(1),FocB1(2),0.2,Colors(3,:),2);
PlotCircle(FocB2(1),FocB2(2),0.2,Colors(3,:),2);
% Apsidenlinie
hp(1)= line([+aB*cos(ups0B) -aB*cos(ups0B)],...
     [+aB*sin(ups0B) -aB*sin(ups0B)],...
     'color',Colors(3,:),'Linewidth',1,'LineStyle',Style(1));
% Start- und Endwert (Randwerte)
hp(2)= line([FocB1(1) x1B],[FocB1(2) y1B],'color',Colors(5,:),...
       'Linewidth',2,'LineStyle',Style(2));
hp(3)= line([FocB1(1) x2B],[FocB1(2) y2B],'color',Colors(5,:),...
       'Linewidth',2,'LineStyle',Style(3));
lgd =legend(hp,'Apsidenlinie','r_1, t_1','r_2, t_2','location','east');
legend box off
set(lgd,'FontSize',12,'FontWeight','normal');
grid on
ylabel('y in 1000 km'); xlabel('x in 1000 km')
axis equal
xlim([-bB bB]*1.2); ylim([-aB aB]*1.2);
ttl=title('Lambert Problem');
set(ttl,'FontSize',14,'FontWeight','normal');
set(gca,'FontSize',16);

%-------------------------------------------------------------------------%
% Ende Programm
% ------------------------------------------------------------------------%

%% Funktionen

% Koordinaten Ellipse im x-y-KOS
function [xe,ye] = xell(a, b, t)
  ye=b.*sin(t); 
  xe=a.*cos(t);
end

% Rotation Ellipse im x-y-KOS um ups0 
function [x1B,y1B]= rotatiere(x1B,y1B,ups0)
    vecin  = [x1B ; y1B; 0]; vecout = mtimes(R_z(-ups0),vecin);
    x1B  = vecout(1); y1B = vecout(2);
end

% Anteil Flugbahn
function xcorr = flugbahn(r1, r2, ecc, a, x, y)
   xcorr = x;
   k=1;
   kend=1;
   while (x(k)-a*ecc)^2 + y(k)^2 <= 0.99*r1^2 
            xcorr(k) = NaN;
            kend = k;
            k = k+1;        
   end
   k = kend-1;
   for k=kend:length(x)
       if (x(k)-a*ecc)^2 + y(k)^2 <= 0.99*r2^2 
            xcorr(k) = NaN;
       end
   end
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
