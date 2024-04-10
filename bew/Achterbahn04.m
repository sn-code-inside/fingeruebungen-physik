% -------------------------------------------------------------------------
% Achterbahn04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet die zeitliche Entwicklung eines Zuges aus einzelnen
% Wagen, die mit Federn verbunden sind, auf einer Gauss-Kurve unter
% Einfluss der Gravitation.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% Parameter
H      = 30;                         % Höhe Berg
B      = 10;                         % Breite Berg
a      = sqrt(log(2))/B;             % Gauss-Parameter 
g      = 9.81;                       % Schwerebeschleunigung
NWagen = 10;                         % Zahl der Wagen
D      = 100000;                     % Federrichtgröße
s0     = 1;                          % Ruhelänge der Federn
m      = 500;                        % Masse eines Wagens
mp     = 70;                         % Masse einer Person

NPoints = 1000;
xw     = linspace(-4/a,4/a,NPoints);
zw     = H*exp(-xw.^2*a^2);
tmax   = 3*sqrt(2*H/g);
tspan  = linspace(0,tmax,NPoints);

% Anfangsbedingungen
Y0  =  zeros(4*NWagen,1);
for i = 1:NWagen
    ix  = (i-1)*4+1; % Index fuer die Variable x von Wagen i
    ixp = (i-1)*4+2; % x'
    iz  = (i-1)*4+3; % z
    izp = (i-1)*4+4; % z'
    Y0(ix)  = -4/a - (i-1)*s0;
    Y0(iz)  = gauss(Y0(ix),a,H);
    hs      = hp(Y0(ix),a,H);
    Y0(ixp) = 1.018*sqrt(2*g*(H-Y0(iz))/(1+hs^2));
    Y0(izp) = hs*Y0(ixp);
end

% Integration der Differentialgleichungen
opt=odeset('AbsTol',1.e-9,'RelTol',1.e-7);           
[t,Y]=ode45(@(t,Y)DGL(t,Y,a,H,g,NWagen,m,D,s0),tspan,Y0,opt); 

%% Abstände und Geschwindigkeiten
figure();
subplot(1, 2, 1);
% Abstand zwischen erstem und zweitem Wagen
diff1(:) = sqrt((Y(:,1)-Y(:,5)).^2+(Y(:,3)-Y(:,7)).^2);
% Abstand zwischen vorleztem und letztem Wagen
diff2(:) = sqrt((Y(:,4*NWagen-3)-Y(:,4*NWagen-7)).^2 ...
    +(Y(:,4*NWagen-1)-Y(:,4*NWagen-5)).^2);
plot(t(:),diff1(:), 'color', Colors(1,:),'LineWidth',2);
hold on
plot(t(:),diff2(:), 'color', Colors(2,:),'LineWidth',2);
%axis([0, tmax, 0, 0.5]);
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
v1(:) = sqrt(Y(:,2).^2+Y(:,4).^2);
% Geschwindigkeit des zweiten Wagens
v2(:) = sqrt(Y(:,4*NWagen-2).^2+Y(:,4*NWagen).^2);
plot(t(:),v1(:), 'color', Colors(3,:),'LineWidth',2);
hold on
plot(t(:),v2(:), 'color', Colors(4,:),'LineWidth',2);
%axis([0, tmax, 0, 0.5]);
ylabel('v in m/s','FontSize',14)
xlabel('t in s','FontSize',14)
grid on
legend('Wagen 1', 'Wagen '+string(NWagen),'location','northwest', ...
    'numcolumns',1);
legend box off
ttl=title('Geschwindigkeit einzelner Wagen über Zeit');
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);


%% Geschwindigkeiten über dem Ort x
w1 = 1; % Nummern der drei Wagen, die betrachtet werden
w2 = 5;
w3 = 10;
% Geschwindigkeiten der drei Wagen
v1(:) = sqrt(Y(:,4*(w1-1)+2).^2+Y(:,4*(w1-1)+4).^2);
v2(:) = sqrt(Y(:,4*(w2-1)+2).^2+Y(:,4*(w2-1)+4).^2);
v3(:) = sqrt(Y(:,4*(w3-1)+2).^2+Y(:,4*(w3-1)+4).^2);

figure();
plot(Y(:,4*(w1-1)+1),v1(:), 'color', Colors(1,:),'LineWidth',2);
hold on
plot(Y(:,4*(w2-1)+1),v2(:), 'color', Colors(2,:),'LineWidth',2);
hold on
plot(Y(:,4*(w3-1)+1),v3(:), 'color', Colors(3,:),'LineWidth',2);
ylabel('v in m/s','FontSize',14)
xlabel('x in m','FontSize',14)
plot(xw,zw)
grid on
legend('Wagen '+string(w1),'Wagen '+string(w2),'Wagen '+string(w3), ...
       'Bahnkurve','location','southwest','numcolumns',1);
legend box off
ttl=title('Geschwindigkeit der einzelen Wagen am Ort x');
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);


%% Zwangskräfte
% Berechne Zwangskräfte der Schienen auf die Wagen
zk1 = Zwangskraefte(w1,Y,a,H,g,NWagen,m,D,s0);
zk2 = Zwangskraefte(w2,Y,a,H,g,NWagen,m,D,s0);
zk3 = Zwangskraefte(w3,Y,a,H,g,NWagen,m,D,s0);
% Berechne Zwangskräfte auf die Personen im Wagen
zkp1 = ZwangskraeftePersonen(w1,Y,a,H,g,mp);
zkp2 = ZwangskraeftePersonen(w2,Y,a,H,g,mp);
zkp3 = ZwangskraeftePersonen(w3,Y,a,H,g,mp);

figure();
subplot(1, 2, 1);
title('Zwangskräfte der Schienen auf die Wagen');
plot(Y(:,4*(w1-1)+1),zk1(:,3)/m/g, 'color', Colors(1,:),'LineWidth',2);
hold on
plot(Y(:,4*(w2-1)+1),zk2(:,3)/m/g, 'color', Colors(2,:),'LineWidth',2);
hold on
plot(Y(:,4*(w3-1)+1),zk3(:,3)/m/g, 'color', Colors(3,:),'LineWidth',2);
ylabel('F/(mg)','FontSize',14)
xlabel('x in m','FontSize',14)
grid on
plot(xw,zw*max(zkp1(:,3)/mp/g/H));
legend('Wagen '+string(w1),'Wagen '+string(w2),'Wagen '+string(w3), ...
       'Bahnkurve','location','southwest','numcolumns',1);
legend box off
ttl=title('Zwangskräfte der Schienen auf die Wagen');
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);

subplot(1, 2, 2);
title('Zwangskräfte auf Personen im Wagen');
plot(Y(:,4*(w1-1)+1),zkp1(:,3)/mp/g, 'color', Colors(1,:),'LineWidth',2);
hold on
plot(Y(:,4*(w2-1)+1),zkp2(:,3)/mp/g, 'color', Colors(2,:),'LineWidth',2);
hold on
plot(Y(:,4*(w3-1)+1),zkp3(:,3)/mp/g, 'color', Colors(3,:),'LineWidth',2);
ylabel('F/(mg)','FontSize',14)
xlabel('x in m','FontSize',14)
grid on
plot(xw,zw*max(zkp1(:,3)/mp/g/H));
legend('Wagen '+string(w1),'Wagen '+string(w2),'Wagen '+string(w3),...
       'Bahnkurve', 'location','southwest','numcolumns',1);
legend box off
ttl=title('Zwangskräfte auf Personen im Wagen');
set(ttl, 'FontSize',12, 'FontWeight' ,'normal');
set(gca,'FontSize',16);
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
function  dY  =  DGL(t,Y,a,H,g,N,m,D,s0)
dY  =  zeros(4*N,1);  %  Es muss ein Spaltenvektor zurückgegeben werden 

% Wagen 1
hs     = hp(Y(1),a,H);
hss    = h2p(Y(1),a,H);
diff   = sqrt((Y(5)-Y(1))^2+(Y(7)-Y(3))^2);
lambda = (m*g+m*hss*Y(2)^2+D*(diff-s0)/diff* ...
    (hs*(Y(5)-Y(1))-(Y(7)-Y(3))))/(1+hs^2);
dY(1)  = Y(2);
dY(2)  = (D*(diff-s0)/diff*(Y(5)-Y(1))-lambda*hs)/m; 
dY(3)  = Y(4);
dY(4)  = (D*(diff-s0)/diff*(Y(7)-Y(3))-m*g+lambda)/m; 


% Wagen 2 bis N-1
for i = 2:N-1
    ix  = (i-1)*4+1; % Index fuer die Variable x von Wagen i
    ixp = (i-1)*4+2; % x'
    iz  = (i-1)*4+3; % z
    izp = (i-1)*4+4; % z'
    hs     = hp(Y(ix),a,H);
    hss    = h2p(Y(ix),a,H);
    diff1  = sqrt((Y(ix)-Y(ix-4))^2+(Y(iz)-Y(iz-4))^2);
    diff2  = sqrt((Y(ix+4)-Y(ix))^2+(Y(iz+4)-Y(iz))^2);
    lambda = (m*g+m*hss*Y(ixp)^2-D*(diff1-s0)/diff1* ...
    (hs*(Y(ix)-Y(ix-4))-(Y(iz)-Y(iz-4)))+D*(diff2-s0)/diff2* ...
    (hs*(Y(ix+4)-Y(ix))-(Y(iz+4)-Y(iz))))/(1+hs^2);
    dY(ix)  = Y(ixp);
    dY(ixp) = (-D*(diff1-s0)/diff1*(Y(ix)-Y(ix-4)) ...
        +D*(diff2-s0)/diff2*(Y(ix+4)-Y(ix))-lambda*hs)/m; 
    dY(iz)  = Y(izp);
    dY(izp) = (-D*(diff1-s0)/diff1*(Y(iz)-Y(iz-4)) ...
        +D*(diff2-s0)/diff2*(Y(iz+4)-Y(iz))-m*g+lambda)/m; 
end

% Wagen N
ix  = (N-1)*4+1; % Index fuer die Variable x von Wagen i
ixp = (N-1)*4+2; % x'
iz  = (N-1)*4+3; % z
izp = (N-1)*4+4; % z'
hs     = hp(Y(ix),a,H);
hss    = h2p(Y(ix),a,H);
diff  = sqrt((Y(ix)-Y(ix-4))^2+(Y(iz)-Y(iz-4))^2);
lambda = (m*g+m*hss*Y(ixp)^2-D*(diff-s0)/diff* ...
    (hs*(Y(ix)-Y(ix-4))-(Y(iz)-Y(iz-4))))/(1+hs^2);
dY(ix)  = Y(ixp);
dY(ixp) = (-D*(diff-s0)/diff*(Y(ix)-Y(ix-4))-lambda*hs)/m; 
dY(iz)  = Y(izp);
dY(izp) = (-D*(diff-s0)/diff*(Y(iz)-Y(iz-4))-m*g+lambda)/m;
end


% Normalkräfte der Schiene auf die Wagen
function  zk = Zwangskraefte(w,Y,a,H,g,NWagen,m,D,s0)

% Fallunterscheidung
if w == 1
    hs     = hp(Y(:,1),a,H);
    hss    = h2p(Y(:,1),a,H);
    diff   = sqrt((Y(:,5)-Y(:,1)).^2+(Y(:,7)-Y(:,3)).^2);
    lambda = (m*g+m.*hss(:).*Y(:,2).^2+D.*(diff(:)-s0)./diff(:).* ...
        (hs(:).*(Y(:,5)-Y(:,1))-(Y(:,7)-Y(:,3))))./(1+hs(:).^2);
    zkx     = -lambda(:).*hs(:); 
    zkz     = lambda(:); 
elseif w == NWagen
    ix  = (NWagen-1)*4+1; % Index fuer die Variable x von Wagen w
    ixp = (NWagen-1)*4+2; % x'
    iz  = (NWagen-1)*4+3; % z
    izp = (NWagen-1)*4+4; % z'
    hs     = hp(Y(:,ix),a,H);
    hss    = h2p(Y(:,ix),a,H);
    diff   = sqrt((Y(:,ix)-Y(:,ix-4)).^2+(Y(:,iz)-Y(:,iz-4)).^2);
    lambda = (m*g+m.*hss(:).*Y(:,ixp).^2-D.*(diff(:)-s0)./diff(:).* ...
    (hs(:).*(Y(:,ix)-Y(:,ix-4))-(Y(:,iz)-Y(:,iz-4))))./(1+hs(:).^2);
    zkx    = -lambda(:).*hs(:); 
    zkz    = lambda(:);
else
    ix  = (w-1)*4+1; % Index fuer die Variable x von Wagen w
    ixp = (w-1)*4+2; % x'
    iz  = (w-1)*4+3; % z
    izp = (w-1)*4+4; % z'
    hs     = hp(Y(:,ix),a,H);
    hss    = h2p(Y(:,ix),a,H);
    diff1  = sqrt((Y(:,ix)-Y(ix-4)).^2+(Y(iz)-Y(iz-4)).^2);
    diff2  = sqrt((Y(:,ix+4)-Y(:,ix)).^2+(Y(:,iz+4)-Y(:,iz)).^2);
    lambda = (m*g+m.*hss(:).*Y(:,ixp).^2-D.*(diff1(:)-s0)./diff1(:).* ...
        (hs(:).*(Y(:,ix)-Y(:,ix-4))-(Y(:,iz)-Y(:,iz-4))) ...
        +D.*(diff2(:)-s0)./diff2.*(hs(:).*(Y(:,ix+4)-Y(:,ix)) ...
        -(Y(:,iz+4)-Y(:,iz))))./(1+hs(:).^2);
    zkx    = -lambda(:).*hs(:); 
    zkz    = lambda(:); 
end

zk(:,1) = zkx(:);
zk(:,2) = zkz(:);
zk(:,3) = sign(zkz(:)).*sqrt(zkx(:).^2+zkz(:).^2);
end


% Normalkräfte, die von den Personen in den Wagen empfunden werden
function  zkp = ZwangskraeftePersonen(w,Y,a,H,g,mp)

ix  = (w-1)*4+1; % Index fuer die Variable x von Wagen w
ixp = (w-1)*4+2; % x'
iz  = (w-1)*4+3; % z
izp = (w-1)*4+4; % z'
hs     = hp(Y(:,ix),a,H);
hss    = h2p(Y(:,ix),a,H);
lambda = (g+hss(:).*Y(:,ixp).^2)./(1+hs(:).^2)*mp;
zkx    = -lambda(:).*hs(:); 
zkz    = lambda(:); 

zkp(:,1) = zkx(:);
zkp(:,2) = zkz(:);
zkp(:,3) = sign(zkz(:)).*sqrt(zkx(:).^2+zkz(:).^2);
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
