% -------------------------------------------------------------------------
% Achterbahn06.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung eines Zuges aus einzelnen
% Wagen, die mit Federn verbunden sind, auf einer Kreisbahn. Dasselbe
% Problem wird für ein kontinuierliches Zug-Modell gelöst.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


% Parameter
R       = 10;                         % Radius des Kreises
g       = 9.81;
NWagen  = 5;                          % Zahl der Wagen
D       = 100000000;                  % Federrichtgröße
s0      = 1.0;                        % Ruhelänge der Federn
m       = 500;                        % Masse eines Wagens
mp      = 70;                         % Masse einer Person
NPoints = 1000;

% Länge des Zuges im kontinuierlichen Modell
dphi = acos(1-0.5*s0^2/R^2); % Winkelabstand zweier Wagen entlang des
                             % Kreisbogens, sodass die gerade gespannten
                             % Federn die Länge s0 haben
L    = (NWagen-1)*R*dphi;    % Länge des Zuges entlang des Kreisbogens

tmax   = 3*sqrt(2*R/g);
tspan  = linspace(0,tmax,NPoints);

%% Lösen der Bewegungsgleichungen im Modell einer Kette aus Wagen
% Anfangsbedingungen
Y0  =  zeros(2*NWagen,1);
for i = 1:NWagen
    iphi      = (i-1)*2+1; % Index fuer die Variable phi von Wagen i
    iphip     = (i-1)*2+2; % d(phi)/dt 
    Y0(iphi)  = (NWagen-i)*dphi; % Winkel, an dem sich Wagen i befindet
    [xw,xwp,zw,zwp] = OrteGeschwindigkeiten(R,Y0(iphi),Y0(iphip));
    Y0(iphip) = sqrt(6*g/R); % alle Wagen sind gleich schnell
end

% Integration der Differentialgleichungen
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);           
[t,Y]=ode45(@(t,Y)DGL(t,Y,R,g,NWagen,m,D,s0),tspan,Y0,opt); 

% Berechnen der Orte und Geschwindigkeiten 
for i = 1:NWagen
    iphi  = (i-1)*2+1; % Index fuer die Variable phi von Wagen i
    iphip = (i-1)*2+2; % d(phi)/dt 
    ix    = (i-1)*4+1;
    ixp   = (i-1)*4+2;
    iz    = (i-1)*4+3;
    izp   = (i-1)*4+4;
    [Pos(:,ix),Pos(:,ixp),Pos(:,iz),Pos(:,izp)] ...
        = OrteGeschwindigkeiten(R,Y(:,iphi),Y(:,iphip));
end

%% Lösen der Bewegungsgleichungen im kontinuierlichen Modell
% Anfangsbedingungen identisch zum ersten Wagen in der Kette
Yk0(1) = Y0(1);
Yk0(2) = Y0(2);

% Integration der Differentialgleichungen
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);           
[tk,Yk]=ode45(@(t,Y)DGLkont(t,Y,L,R,g),tspan,Yk0,opt); 

% Berechnen der Orte und Geschwindigkeiten 
[Posk(:,1),Posk(:,2),Posk(:,3),Posk(:,4)] ...
        = OrteGeschwindigkeiten(R,Y(:,1),Y(:,2));


%% Abstände und Geschwindigkeiten
figure();
subplot(1, 2, 1);
% Abstand zwischen erstem und zweitem Wagen
diff1(:) = sqrt((Pos(:,1)-Pos(:,5)).^2+(Pos(:,3)-Pos(:,7)).^2);
% Abstand zwischen vorleztem und letztem Wagen
diff2(:) = sqrt((Pos(:,4*NWagen-3)-Pos(:,4*NWagen-7)).^2 ...
    +(Pos(:,4*NWagen-1)-Pos(:,4*NWagen-5)).^2);
plot(t(:),diff1(:), 'color', Colors(1,:),'LineWidth',3);
hold on
plot(t(:),diff2(:), 'color', Colors(2,:),'LineWidth',2);
axis([0, tmax, 0.8, 1.2]);
ylabel('Abstand in m','FontSize',14)
xlabel('t in s','FontSize',14)
grid on
legend('Wagen 1-Wagen 2', 'Wagen '+string(NWagen-1)+'-Wagen '+ ...
    string(NWagen),'location','southwest','numcolumns',1);
legend box off
ttl=title('Abstand einzelner Wagen über Zeit');
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);

% Geschwindigkeiten
subplot(1, 2, 2);
% Geschwindigkeit des ersten Wagens
v1(:) = sqrt(Pos(:,2).^2+Pos(:,4).^2);
% Geschwindigkeit des zweiten Wagens
v2(:) = sqrt(Pos(:,4*NWagen-2).^2+Pos(:,4*NWagen).^2);
plot(t(:),R*Yk(:,2), 'color', Colors(6,:),'LineWidth',3);
hold on;
plot(t(:),v1(:), 'color', Colors(3,:),'LineWidth',2);
hold on
plot(t(:),v2(:), 'color', Colors(4,:),'LineWidth',2);
ylabel('v in m/s','FontSize',14)
xlabel('t in s','FontSize',14)
grid on
legend('kontinuierliches Modell', 'Wagen 1', 'Wagen '+string(NWagen), ...
    'location','northwest','numcolumns',1);
legend box off
ttl=title('Geschwindigkeit einzelner Wagen über Zeit');
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktionen

% Orte und Geschwindigkeiten im Ortsraum
function [x,vx,z,vz] = OrteGeschwindigkeiten(R,phi,phip)
    x  = R.*sin(phi);
    vx = phip.*R.*cos(phi);
    z  = R.*(1-cos(phi));
    vz = phip.*R.*sin(phi);
 end


% Differentialgleichungssystem
function  dY  =  DGL(t,Y,R,g,N,m,D,s0)
dY  =  zeros(2*N,1);  %  Es muss ein Spaltenvektor zurückgegeben werden 

% Wagen 1
[x1,xp1,xpp1,z1,zp1,zpp1] = Position(R,Y(1));
[x2,xp2,xpp2,z2,zp2,zpp2] = Position(R,Y(3));
diff   = sqrt((x2-x1)^2+(z2-z1)^2);
dY(1)  = Y(2);
dY(2)  = (-(xp1*xpp1+zp1*zpp1)*Y(2)^2 ...
    +D/m*(diff-s0)/diff*(x2-x1)*xp1 ...
    +D/m*(diff-s0)/diff*(z2-z1)*zp1 - g*zp1)/(xp1^2+zp1^2);

% Wagen 2 bis N-1
for i = 2:N-1
    [xi,xpi,xppi,zi,zpi,zppi] = Position(R,Y((i-1)*2+1)); % i
    [xip1,xpip1,xppip1,zip1,zpip1,zppip1] = Position(R,Y(i*2+1));     %i+1
    [xim1,xpim1,xppim1,zim1,zpim1,zppim1] = Position(R,Y((i-2)*2+1)); % i-1
    iphi  = (i-1)*2+1; % Index fuer die Variable phi von Wagen i
    iphip = (i-1)*2+2; % d(phi)/dt
    diff1     = sqrt((xi-xim1)^2+(zi-zim1)^2);
    diff2     = sqrt((xip1-xi)^2+(zip1-zi)^2);
    dY(iphi)  = Y(iphip);
    dY(iphip) = (-(xpi*xppi+zpi*zppi)*Y(iphip)^2 ...
        -D/m*(diff1-s0)/diff1*(xi-xim1)*xpi ...
        +D/m*(diff2-s0)/diff2*(xip1-xi)*xpi ... 
        -D/m*(diff1-s0)/diff1*(zi-zim1)*zpi ...
        +D/m*(diff2-s0)/diff2*(zip1-zi)*zpi - g*zpi)/(xpi^2+zpi^2);
end

% Wagen N
[xN,xpN,xppN,zN,zpN,zppN] = Position(R,Y((N-1)*2+1));
[xNm1,xpNm1,xppNm1,zNm1,zpNm1,zppNm1] = Position(R,Y((N-2)*2+1));
diff   = sqrt((xN-xNm1)^2+(zN-zNm1)^2);
dY((N-1)*2+1)  = Y((N-1)*2+2);
dY((N-1)*2+2)  = (-(xpN*xppN+zpN*zppN)*Y((N-1)*2+2)^2 ...
    -D/m*(diff-s0)/diff*(xN-xNm1)*xpN ...
    -D/m*(diff-s0)/diff*(zN-zNm1)*zpN - g*zpN)/(xpN^2+zpN^2);
end


% Differentialgleichungssystem für das kontinuierliche Modell
function  dY  =  DGLkont(t,Y,L,R,g)
dY  =  zeros(2,1);  %  Es  muss  ein  Spaltenvektor  zurückgegeben  werden
dY(1) = Y(2);
dY(2) = g/L*(cos(Y(1))-cos(Y(1)-L/R));
end

% Positionen im Ortsraum
function [x,xp,xpp,z,zp,zpp] = Position(R,phi)
    x   = R*sin(phi);
    xp  = R*cos(phi);
    xpp = -R*sin(phi);
    z   = R*(1-cos(phi));
    zp  = R*sin(phi);
    zpp = R*cos(phi);
end

% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------
