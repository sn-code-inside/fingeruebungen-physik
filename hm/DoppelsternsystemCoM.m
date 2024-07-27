%%-------------------------------------------------------------------------
% DoppelsternsystemCoM.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Himmelsmechanik" aus
% "FingerÃ¼bungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% 
% Bewegung eines Doppelsternsystems mit Massenverhältnis 2.5:1
% Vergleich numerischer und analytischer Rechnung im Schwerpunktsystem
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

day  = 24 * 60 * 60;    % Tag in s
year = 365.25 * day;    % Jahr in s
AE   = 1.495978707e11;  % AE in m
G    = 6.6743e-11;      % G in m³ / (kg * s²)
m0   = 2e30;            % Massen in kg

% Simulationszeit und Zeitschrittweite [in s].
tend  = 20 * year;
dt    = 0.05 * day;


Y0    = 1.0 * AE;     % Abstand des Körpers vom Schwerpunkt.
V0    = 25e3;       % Geschwindigkeit  (Betrag m/s).

%% Anfangspositionen und -geschwindigkeiten der Körper
NK = 2;

vecr01 = [0;-1]*AE;
vecr02 = [0;1]*AE;
vecv01 = [-0.2;0]*V0;
vecv02 = [0.8;0]*V0;
m = [1 2/5]*m0;
q(1) = 3*m(1)/m(1);
q(2) = 3*m(2)/m(1);

%% Schwerpunktsystem

vecrS  = (m(1)*vecr01 + m(2)*vecr02)/(m(1)+m(2));
vecr01 = vecr01 - vecrS;
vecr02 = vecr02 - vecrS;
vecvS  = (m(1)*vecv01 + m(2)*vecv02)/(m(1)+m(2));
vecv01 = vecv01 - vecvS;
vecv02 = vecv02 - vecvS;


%% ode45

AB = [vecr01, vecr02, vecv01, vecv02];
P1.G   = G;
P1.m1  = m(1);
P1.m2  = m(2);
P1.NK  = NK;
opts = odeset('AbsTol',1.e-13,'RelTol',1.e-12);
[tA,YA]=ode45(@(t,YA, P1)DGL_3BodySystem(t,YA,P1),[0 tend],AB,opts,P1);
tA = tA/86400;


%% Graphik 

figure('name','Trajektorien')
hold on
axmin = -2;
axmax = 2;
aymin = -2;
aymax = +2;
axis([axmin axmax aymin aymax])

plot(vecr01(1)/AE,vecr01(2)/AE,'o','markersize',...
        2*q(1),'linewidth',2*q(1),'color',Colors(2,:));
plot(vecr02(1)/AE,vecr02(2)/AE,'o','markersize',...
        2*q(2),'linewidth',2*q(2),'color',Colors(3,:));
axis equal
line([vecr01(1) vecr02(1)]/AE,[vecr01(2) vecr02(2)]/AE,...
        'color',Colors(4,:))
grid on

for k=1:NK
    hp(k) =plot(YA(:,2*k-1)/AE,YA(:,2*k)/AE,'Color',Colors(k+1,:),...
    'linewidth',2);
end
hp(3)=plot(0,0,'+k','markersize',10,'linewidth',2);

% Graphik Simulation der Bewegung
h = text(-2.0,1.5,strcat(num2str(0,'%6.0f'),' Tage'));
for n=1:100:length(tA) 
for k=1:NK
    set(h,'Visible','off');
    plot(YA(n,2*k-1)/AE,YA(n,2*k)/AE,'o','markersize',2*q(k),...
        'linewidth',q(k),'color',Colors(k+1,:))
    h = text(-2,1.5,strcat(num2str(tA(n),'%6.0f'),' Tage'));
    set(h,'Visible','on', 'FontSize',12);
    pause(0.05)
end
end
legend(hp(1:3),'m_1', 'm_2', 'Schwerpunkt','location', 'northeast')
legend box off
xlabel('x in AE')
ylabel('y in AE')
set(gca,'FontSize',14);


%% Energie, Drehimpuls und Schwerpunkt
x(:,1)  = YA(:,1);
y(:,1)  = YA(:,2);
x(:,2)  = YA(:,3);
y(:,2)  = YA(:,4);
vx(:,1) = YA(:,5);
vy(:,1) = YA(:,6);
vx(:,2) = YA(:,7);
vy(:,2) = YA(:,8);
dr12(:) = sqrt((x(:,1)- x(:,2)).^2 + (y(:,1)- y(:,2)).^2);

XS(:) = (m(1)*x(:,1)+m(2)*x(:,2))/(m(1)+m(2));
YS(:) = (m(1)*y(:,1)+m(2)*y(:,2))/(m(1)+m(2));

% Berechne die verschiedenen Energiebeiträge.
Ekin(:,1) = vx(:,1).^2+ vy(:,1).^2;
Ekin(:,2) = vx(:,2).^2+ vy(:,2).^2;
Ekinges(:) = (m(1)*Ekin(:,1)+m(2)*Ekin(:,2))/2;
Epot(:) = -G * m(1)*m(2) *(1./dr12(:));
E0 = Ekinges(1)+Epot(1);

TJ =1e12;
figure('Name','Energiebeiträge')
yyaxis left
hp2(1)=plot(tA, Ekinges/TJ,'linewidth',2,'linestyle',Style(2));
hold on
hp2(2)=plot(tA, Epot/TJ,'linewidth',2,'linestyle',Style(3));
hp2(3)=plot(tA, (Ekinges+Epot)/TJ,'linewidth',2,'linestyle',Style(1));
ylabel('E in TJ')
yyaxis right
hp2(4)=plot(tA, 1-(Ekinges+Epot)/E0,'linewidth',2,'linestyle',Style(1));
ylim([0 1.e-11])
grid on
legend(hp2,'T_{kin}', 'U_{pot}', 'E_{ges}','1-E(t)/E_0','location', 'best')
legend box off
ylabel('1-E(t)/E_0')
xlabel('t in Tagen')
set(gca,'FontSize',14);

L1(:,:) = m(1)*cross([x(:,1), y(:,1), y(:,1)*0],[vx(:,1), vy(:,1), vy(:,1)*0]);
L2(:,:) = m(2)*cross([x(:,2), y(:,2), y(:,2)*0],[vx(:,2), vy(:,2), vy(:,2)*0]);
LG(:) = L1(:,3) + L2(:,3);

figure('Name','Drehimpulsbeiträge und Schwerpunkt')
yyaxis left
hp2(1)=plot(tA, 1-LG/LG(1),'linewidth',2, 'color', Colors(3,:));
ylabel('1-L(t)/L_0')
yyaxis right
hp2(2)=plot(tA,XS,'linewidth',2,'color',Colors(4,:),'linestyle',Style(1));
hold on
hp2(3)=plot(tA,YS,'linewidth',2,'color',Colors(4,:),'linestyle',Style(2));
ylim([-2e-2 1e-2])

grid on
legend(hp2(1:3),'1-L(t)/L_0','x_S', 'y_S', 'location', 'southwest')
legend box off
ylabel('Schwerpunktkoordinaten in m')
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


