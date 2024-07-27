%%-------------------------------------------------------------------------
% Doppelsternsystem.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Himmelsmechanik" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Bewegung eines Doppelsternsystems mit Massenverhältnis 2.5:1
% Vergleich numerischer und analytischer Rechnung
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
day  = 86400; % Tag in s
year = 365.25 * day; % Jahr in s
AE   = 1.495978707e11;  % AE in km
G    = 6.6743e-11;      % G in m³ / (kg * s²)
m0   = 2e30;            % Massen in kg

% Simulationszeit und Zeitschrittweite [in s].
tend  = 5 * year;
dt    = 0.1 * day;


Y0    = 1 * AE;     % Abstand des Körpers vom Schwerpunkt.
V0    = 25e3;       % Geschwindigkeit  (Betrag m/s).

%%  Anfangspositionen und -geschwindigkeiten der Körper
NK = 2;
vecr01 = [0;-1]*AE;
vecr02 = [0;1]*AE;
vecv01 = [-0.2;0]*V0;
vecv02 = [0.8;0]*V0;
m = [1 2/5]*m0;
q(1) = 3*m(1)/m(1);
q(2) = 3*m(2)/m(1);
AB = [vecr01, vecr02, vecv01, vecv02];
vecrS0  = (m(1)*vecr01 + m(2)*vecr02)/(m(1)+m(2));

P1.G   = G;
P1.m1  = m(1);
P1.m2  = m(2);
P1.NK  = NK;

%% ode45 

opts = odeset('AbsTol',1.e-8,'RelTol',1.e-9);
[tA,YA]=ode45(@(t,YA, P1)DGL_3BodySystem(t,YA,P1),[0 tend],AB,opts,P1);
tA = tA/day;

%% Graphik Trajektorien
figure('name','Trajektorien')
hold on
axmin = -1.0;
axmax = 4;
aymin = -2.5;
aymax = +2.5;
hold on
axis([axmin axmax aymin aymax])
axis equal

hp(1)=plot(vecr01(1)/AE,vecr01(2)/AE,'o','markersize',...
        2*q(1),'linewidth',2*q(1),'color',Colors(2,:));
hp(2)=plot(vecr02(1)/AE,vecr02(2)/AE,'o','markersize',...
        2*q(2),'linewidth',2*q(2),'color',Colors(3,:));
hp(3)=plot(vecrS0(1)/AE,vecrS0(2)/AE,'+k','markersize',...
        12,'linewidth',2);
axis equal
line([vecr01(1) vecr02(1)]/AE,[vecr01(2) vecr02(2)]/AE,...
        'colo',Colors(4,:))
grid on

for k=1:NK
    plot(YA(:,2*k-1)/AE,YA(:,2*k)/AE,'Color',Colors(k+1,:))
end

%% Schwerpunktsystem

vecrS(:,:)  = (m(1)*[YA(:,1) YA(:,2)]+ m(2)*[YA(:,3) YA(:,4)])/(m(1)+m(2));

h = text(0,2,strcat(num2str(0,'%6.0f'),' Tage'));
for n=1:10:length(tA) 
for k=1:NK
    set(h,'Visible','off');
    plot(YA(n,2*k-1)/AE,YA(n,2*k)/AE,'o','markersize',q(k),...
        'linewidth',q(k),'color',Colors(k+1,:))
    plot(YA(n,2*k-1)/AE,YA(n,2*k)/AE,'o','markersize',q(k),...
        'linewidth',q(k),'color',Colors(k+1,:))
    line([vecrS(1,1) vecrS(n,1)]/AE, [vecrS(1,2) vecrS(n,2)]/AE,...
          'color','k','linewidth',1)
    h = text(0,2.0,strcat(num2str(tA(n),'%6.0f'),' Tage'));
    set(h,'Visible','on', 'FontSize',12);
    pause(0.01)
end
end
plot(vecrS(1,1)/AE,vecrS(1,2)/AE,'+k','markersize',12,'linewidth',2)
plot(vecrS(end,1)/AE,vecrS(end,2)/AE,'+k','markersize',12,'linewidth',2)
legend(hp,'m_1', 'm_2','Schwerpunkt','location', 'northeast')
legend box off
axis([axmin axmax aymin aymax])
axis equal
set(gca,'FontSize',14);
xlabel('x in AE')
ylabel('y in AE')


x(:,1)  = YA(:,1);
y(:,1)  = YA(:,2);
x(:,2)  = YA(:,3);
y(:,2)  = YA(:,4);
vx(:,1) = YA(:,5);
vy(:,1) = YA(:,6);
vx(:,2) = YA(:,7);
vy(:,2) = YA(:,8);
dr12(:) = sqrt((x(:,1)- x(:,2)).^2 + (y(:,1)- y(:,2)).^2);

% Berechne die verschiedenen Energiebeiträge.
Ekin(:,1) = vx(:,1).^2+ vy(:,1).^2;
Ekin(:,2) = vx(:,2).^2+ vy(:,2).^2;
Ekinges(:) = (m(1)*Ekin(:,1)+m(2)*Ekin(:,2))/2;
 
Epot(:) = -G * m(1)*m(2) *(1./dr12(:));
TJ =1e12;
figure('Name','Energiebeiträge')
hp2(1)=plot(tA, Ekinges/TJ,'linewidth',2);
hold on
hp2(2)=plot(tA, Epot/TJ,'linewidth',2);
hp2(3)=plot(tA, (Ekinges+Epot)/TJ,'linewidth',2);
grid on
legend(hp2,'T_{kin}', 'U_{pot}', 'E_{ges}', 'location', 'best')
legend box off
ylabel('E in TJ')
xlabel('t in Tagen')
set(gca,'FontSize',14);

Zeros(:) = transpose((x(:,1)-x(:,1)));
L1(:,:) = m(1)*cross([x(:,1), y(:,1), y(:,1)*0],[vx(:,1), vy(:,1), vy(:,1)*0]);
L2(:,:) = m(2)*cross([x(:,2), y(:,2), y(:,2)*0],[vx(:,2), vy(:,2), vy(:,2)*0]);
LG(:) = L1(:,3) + L2(:,3);

figure('Name','Drehimpulsbeiträge')
hp2(1)=plot(tA, L1(:,3),'linewidth',2, 'color', Colors(2,:));
hold on
hp2(2)=plot(tA, L2(:,3),'linewidth',2, 'color', Colors(3,:));
hp2(3)=plot(tA, LG,'linewidth',2, 'color', Colors(4,:));
grid on
legend(hp2,'L_1', 'L_2', 'L_{ges}', 'location', 'southwest')
legend box off
ylabel('L in kg m²/s')
xlabel('t in Tagen')
set(gca,'FontSize',14);


%%  Funktion
% DGL
function dY = DGL_3BodySystem(~, Y, P1)
    G   = P1.G;
    m1  = P1.m1;
    m2  = P1.m2;
    NK  = P1.NK;
    for k=1:NK
      vecr(:,k) = [Y(2*k-1) Y(2*k)];
      vecv(:,k) = [Y(2*k+3) Y(2*k+4)];
      x(k)      = vecr(1,k);
      y(k)      = vecr(2,k);
    end  
    dr12 = norm(vecr(:,1)-vecr(:,2));
    dY = [Y(5);...
          Y(6);...
          Y(7);...
          Y(8);...
          -G*m2*(x(1)-x(2))/dr12^3;...
          -G*m2*(y(1)-y(2))/dr12^3;...
          -G*m1*(x(2)-x(1))/dr12^3;...
          -G*m1*(y(2)-y(1))/dr12^3       
         ];
end

%% Ende Programm


