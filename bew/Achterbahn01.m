% -------------------------------------------------------------------------
% Achterbahn01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung eines Massenpunktes auf
% einer Gauss-Kurve unter Einfluss der Gravitation.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Anfangsbedingungen, Parameter
H  = 30;                             % Höhe Berg
B  = 10;                             % Breite Berg
a  = sqrt(log(2))/B;                 % Gauss-Parameter 
g  = 9.81;                           % Schwerebeschleunigung

NPoints = 1000;
xw     = linspace(-4/a,4/a,NPoints);
tspan  = linspace(0,3*sqrt(2*H/g),NPoints);

x0   = -4/a;
hsx0 = -2*H*a^2*x0*exp(-x0^2*a^2);
dx0  = 1.018*sqrt(2*g*H/1+hsx0^2);
z0   = H*exp(-x0^2*a^2);
dz0  = hsx0*dx0;

AB = [x0, dx0, z0, dz0];
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);           
[t,Y]=ode45(@(t,Y)DGL(t,Y,a,H,g),tspan,AB,opt); 
x  = Y(:,1);
vx = Y(:,2);
z  = Y(:,3);
vz = Y(:,4);
T  = 0.5*(vx.^2+vz.^2);
U  = g*z;

% Ableitungen
hs  = hp(x,a, H);
h2s = h2p(x,a,H);
fac = 1./(1+hs.^2);

% Beschleunigungen
ax = -(h2s.*vx.^2+g).*hs.*fac;
az = h2s.*vx.^2 + hs.*ax ;

% Normierung für graphische Darstellung
agx = ax*H/10/g;
agz = az*H/10/g;
Nz  = (az+g);
Ngz = Nz*H/10/g;
Ngx = agx;
N   = sqrt(Nz.^2+ax.^2)/g;
v   = sqrt(vz.^2+vx.^2);
theta1 = acosd((ax.*vx+Nz.*vz)./(N.^2+v.^2));
theta2 = atand(hs);

%% Bewegung auf Bahnkurve und Kräfte (Simulation)
figure()
hold on;
h=title('Achterbahn');
set(h,'FontSize',14,'FontWeight','normal'); 
axis([x0, -x0, -0.2*H, 1.2*H]);
ylabel('Höhe h in m','FontSize',14)
xlabel('x in m','FontSize',14)
set(gca,'FontSize',16);
grid on
plot(x,z,'color', Colors(7,:),'LineWidth',2);
for k = 1:NPoints/25:length(t)
        p1 = plot(x(k),z(k),'s','color', Colors(4,:),'LineWidth',2);
      %Geschwindigkeit
      p2 = line([x(k),x(k)+vx(k)/2],[z(k),z(k)+vz(k)/2],...
                    'color', Colors(3,:),'LineWidth',2);
%     %Beschleunigung
%         p2 = line([x(k),x(k)+agx(k)],[z(k),z(k)+agz(k)],...
%                         'color', Colors(4,:),'LineWidth',2);
%     %Normalkraft
%         p3 = line([x(k),x(k)+agx(k)],[z(k),z(k)+Ngz(k)],...
%                         'color', Colors(9,:),'LineWidth',2);
    pause(0.1);
    p1. Visible = 'off';
%     p2. Visible = 'off';
%     p3. Visible = 'off';
    grid on
end
legend(p2,'Geschwindigkeit',...
                'location','northeast','numcolumns',1);
% legend(p3, 'Normalbeschleunigung',...
%                'location','northeast','numcolumns',1);
legend box off
set(gca,'FontSize',16);

%% Energie
figure()  % Energie
Eg = (T(1)+U(1))/100;
plot(x,T/Eg, 'color', Colors(2,:),'LineWidth',2);
hold on
plot(x,U/Eg,'color', Colors(3,:),'LineWidth',2);
plot(x,(T+U)/Eg,'color', Colors(4,:),'LineWidth',2);
axis([x0, -x0, -20, 120]);
ylabel('Energie in % ','FontSize',14)
xlabel('x in m','FontSize',14)
grid on
legend('kin. Energie', 'pot. Energie','Gesamtenergie', ...
               'location','east','numcolumns',1);
legend box off
h=title('Achterbahn - Energie');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

%% Beschleunigungen
figure()
h=title('Achterbahn - Beschleunigung ');
set(h,'FontSize',14,'FontWeight','normal'); 
% yyaxis left
hold on;
% ylabel('Geschwindigkeit in m/s','FontSize',14)
% plot(t,vx,'color', Colors(3,:),'LineWidth',2);
% plot(t,vz,'color', Colors(4,:),'LineWidth',2);
ylabel('Normalkraft in m x g','FontSize',14)
p(1)=plot(x,N,'color', Colors(3,:),'LineWidth',2);
axis([x0, -x0, -2, 8]);
xlabel('x in m','FontSize',14)
grid on
yyaxis right
% hold on
% plot(t,x,'color', Colors(7,:),'LineWidth',2);
axis([x0, -x0, -10, 40]);
p(2)=plot(x,z,'color', Colors(7,:),'LineWidth',2);
grid on
ylabel('Höhe in m','FontSize',14)
set(gca,'FontSize',16);
legend(p,'Normalkraft', 'z(t)',  ...
               'location','south','numcolumns',2);
legend box off
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

%%
% Funktionen

function zg = gauss(x,a,H)
    zg = H*exp(-a^2*x.^2);
end

function hprime = hp(x,a,H)
    hprime = -2*a^2*x.*gauss(x,a,H);
end

function h2prime = h2p(x,a,H)
    h2prime = 2*a^2*(2*a^2*x.^2-1).*gauss(x,a,H);
end


% Differentialgleichungssystem
function  dY  =  DGL(t,Y,a,H,g)
dY  =  zeros(4,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden 
gaussk = H*exp(-Y(1).^2*a^2);
hs     = -2*a^2*Y(1).*gaussk;
hss    =  2*a^2*(2*a^2*Y(1).^2-1).*gaussk;
dY(1)  =  Y(2);
dY(2)  = -(hss*Y(2)^2+g)*hs/(1+hs^2); 
dY(3)  =  Y(4);
dY(4)  = (hss*Y(2)^2-g*hs^2)/(1+hs^2); 
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------

