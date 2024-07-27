%%-------------------------------------------------------------------------
% DreiKoerproblem02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Himmelsmechanik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Bewegung von drei Körpern, mit Anfangsbedingung auf 
% einem gleichseitigen Dreieck. 
% 
% -------------------------------------------------------------------------

%% Initialisierung
clc
clear 
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors=GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter

RE  = 6371;         % in km

day  = 84600;           % Tag in s
year = 365.25 * day;    % Jahr in s
AE   = 1.495978707e11;  % AE in km
G    = 6.6743e-11;      % G in m^3 / (kg s^3)
m0   = 2e30;            % Massen in kg

% Simulationszeit und Zeitschrittweite [in s].
tend  = 12* year;
dt    = 1 * day;

d = 1 * AE;                      % Abstand des Körpers vom Schwerpunkt.
v0 = sqrt(G * m0 / d / sqrt(3)); % Geschwindigkeit von m0 (Betrag).

%% Anfangspositionen und -geschwindigkeitender Körper
NK = 3;
for k=1:NK
  Alpha(k) = (k-1)*360/NK+90;
end
epsilon = 0.0;  %kleine Störung

vecr0 = d .* [cosd(Alpha); sind(Alpha)];
vecv0 = v0.* [-sind(Alpha); cosd(Alpha)];
vecr0(:,1)= vecr0(:,1) + epsilon;  %kleine Störung

for k =1:NK
   AB(2*k-1)   = vecr0(1,k);
   AB(2*k)     = vecr0(2,k);
   AB(2*k+5)   = vecv0(1,k);
   AB(2*k+6)   = vecv0(2,k);
end

%% Berechnung Trajektorien über ode45
P1.G   = G;
P1.m0  = m0;
P1.NK  = NK;
opts = odeset('AbsTol',1.e-14,'RelTol',1.e-13);
[tA,YA]=ode45(@(t,YA, P1)DGL_3BodySystem(t,YA,P1),[0 tend],AB,opts,P1);

%% Graphische Ausgabe

figure('name','Trajektorien')
hold on
for k=1:NK
  hp(k) = plot(vecr0(1,k)/AE,vecr0(2,k)/AE,'o','markersize',8,'linewidth',8,...
       'color',Colors(k,:));
end
axis equal
al =3;
axis([-al al -al al])
hold on
for k=1:length(Alpha)-1
    line([vecr0(1,1) vecr0(1,k+1)]/AE,[vecr0(2,1) vecr0(2,k+1)]/AE,...
        'color',Colors(4,:))
end
line([vecr0(1,3) vecr0(1,2)]/AE,[vecr0(2,3) vecr0(2,2)]/AE,...
    'color',Colors(4,:))
grid on
for k=1:NK
    plot(YA(:,2*k-1)/AE,YA(:,2*k)/AE,'Color',Colors(k,:))
end
for n=1:100:length(tA)
for k=1:NK
    % plot(YA(n,2*k-1)/AE,YA(n,2*k)/AE,'o','markersize',4,'linewidth',4,...
       % 'color',Colors(k,:));
    % pause(0.01)
end
end
legend(hp,'m_1', 'm_2', 'm_3' ,'location', 'northeast')
legend box off
set(gca,'FontSize',14);
xlabel('x in AE')
ylabel('y in AE')


%% Berechnung und Darstellung der Energiebeiträge.

x(:,1) = YA(:,1);
y(:,1) = YA(:,2);
x(:,2) = YA(:,3);
y(:,2) = YA(:,4);
x(:,3) = YA(:,5);
y(:,3) = YA(:,6);
vx(:,1) = YA(:,7);
vy(:,1) = YA(:,8);
vx(:,2) = YA(:,9);
vy(:,2) = YA(:,10);
vx(:,3) = YA(:,11);
vy(:,3) = YA(:,12);


dr12(:) = sqrt((x(:,1)- x(:,2)).^2 + (y(:,1)- y(:,2)).^2);
dr13(:) = sqrt((x(:,1)- x(:,3)).^2 + (y(:,1)- y(:,3)).^2);
dr23(:) = sqrt((x(:,2)- x(:,3)).^2 + (y(:,2)- y(:,3)).^2);

Ekin(:,1) = vx(:,1).^2+ vy(:,1).^2;
Ekin(:,2) = vx(:,2).^2+ vy(:,2).^2;
Ekin(:,3) = vx(:,3).^2+ vy(:,3).^2;
Ekinges(:) = m0*(Ekin(:,1)+Ekin(:,2)+Ekin(:,3))/2;
 
Epot(:) = -G * m0^2 *(1./dr12(:) +1./dr13(:) + 1./dr23(:));

TJ = 1e12;
tA = tA/86400;
figure('name','Energieverhältnisse')
hp2(1) = plot(tA, Ekinges/TJ,'linewidth',2);
hold on
hp2(2) = plot(tA, Epot/TJ,'linewidth',2);
hp2(3) = plot(tA, (Ekinges+Epot)/TJ,'linewidth',2);
grid on
legend(hp2,'T_{kin}', 'U_{pot}', 'E_{ges}', 'location', 'northwest')
legend box off
ylabel('E in TJ')
xlabel('t in Tagen')
set(gca,'FontSize',14);


%% DGL -------------------------------------------------------------------

function dY = DGL_3BodySystem(~, Y, P1)
    G   = P1.G;
    m0  = P1.m0;
    GM2 = G*m0;
    NK  = P1.NK;
    for k=1:NK
      vecr(:,k) = [Y(2*k-1) Y(2*k)];
      vecv(:,k) = [Y(2*k+5) Y(2*k+6)];
      x(k)      = vecr(1,k);
      y(k)      = vecr(2,k);
    end  
    dr12 = norm(vecr(:,1)-vecr(:,2));
    dr13 = norm(vecr(:,1)-vecr(:,3));
    dr23 = norm(vecr(:,2)-vecr(:,3));
    dY =  [Y(7);...
          Y(8);...
          Y(9);...
          Y(10);...
          Y(11);...
          Y(12);...
          -GM2*(x(1)-x(2))/dr12^3-GM2*(x(1)-x(3))/dr13^3;...
          -GM2*(y(1)-y(2))/dr12^3-GM2*(y(1)-y(3))/dr13^3;...
          -GM2*(x(2)-x(1))/dr12^3-GM2*(x(2)-x(3))/dr23^3;...
          -GM2*(y(2)-y(1))/dr12^3-GM2*(y(2)-y(3))/dr23^3;...        
          -GM2*(x(3)-x(1))/dr13^3-GM2*(x(3)-x(2))/dr23^3;...        
          -GM2*(y(3)-y(1))/dr13^3-GM2*(y(3)-y(2))/dr23^3       
         ];
end


%% Ende Programm


