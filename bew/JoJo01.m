% -------------------------------------------------------------------------
% JoJo01.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet Dynamik des Jo-Jos ohne Handbewegung
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerische Lösung

% Parameter

L       =  1.000;           % Länge Stab in  m
Ri      =  0.009;           % Innenradius in m
Ra      =  0.045;           % Aussenradius in m
rmax    =  0.013;           % maximaler Radius in m 
g       =  9.81;            % g in m/s^2 
m       =  0.040;           % Masse in kg 
x0      =  0;               % Anfangsposition
dx0     =  0.0;             % Anfangsgeschwindigkeit in m/s
d       =  pi*(rmax^2-Ri^2)/L;    % effektive Dicke
muF     =  0.01/L;           % Masse pro Fadenlänge
tmax1   =  5;
zspan   =  linspace(0.0,1.0,100); % z-Bereich


fprintf('\n ');
fprintf('\n L   = %8.2f m',  L);
fprintf('\n m   = %8.2f kg', m);
fprintf('\n Ri  = %8.2f mm', Ri*1000);
fprintf('\n d   = %8.2f mm', d*1000);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB1=[x0,dx0]; % AB für ode45

%%
% Dynamikberechnung für verschiedene Jo-Jos
Para_Str = string(4);
d_auswahl = [0.0001, d, 2*d, 4*d];

P1.g  = g;
P1.Ra = Ra;
P1.Ri = Ri;
P1.L  = L;

omegaB   =  2*sqrt(g*L)/Ra;
tmax2    =  2*pi/omegaB;
tspan1   =  linspace(0.0,tmax1,200);
tspan2   =  linspace(0.0,tmax2,200);

options1 = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent1);

z     = NaN(4,200);
z2    = NaN(4,200);
zr    = NaN(4,200);
phi   = NaN(4,200);
v2    = NaN(4,200);
td    = NaN(4,200);
td2   = NaN(4,200);
FT2   = NaN(4,200);

for k=1:4
    d = d_auswahl(k);
    P1.d = d;
    Para_Str(k,:) = strcat(strcat(' d = ', num2str(d*1000,'%3.1f')),' mm');
    rx  = sqrt(Ri^2+d*(L-zspan)/pi);
    J0  = 0.5*m*Ra^2;
    Jx  = J0 + 0.*muF*rx.^2;
    v(k,:)   = sqrt(2*g*zspan)./sqrt(1+Jx/m./rx.^2);
    a(k,:)   = 2*g*rx.^2./(Ra^2+2*rx.^2);
    FT(k,:)  = m*g*Ra^2./(Ra^2+2*rx.^2);
    %Abrollbereich
    [t, Y, TE, YE, IE] = ode45(@dgl_JoJo1,tspan1, AB1,options1,P1); 
    for k1=1:length(t) 
        td(k,k1) = t(k1); 
        z(k,k1)  = Y(k1,1);
        vend(k) = YE(2);
        tfall(k)= TE; 
    end
    for k1=1:length(z(k,:))
        zr(k,k1) = z(k,length(z(k,:))+1-k1);      
    end
    %Umkehrbereich
    AB2=[+pi/2,-vend(k)/Ri]; % AB für ode45
    options2 = odeset('AbsTol',1.e-7,'RelTol',1.e-5,'events',@MyEvent2);
    [t, Y, TE, YE, IE] = ode45(@dgl_JoJo2,tspan2, AB2,options2,P1); 
    for k1=1:length(t) 
        td2(k,k1) = t(k1); 
        phi(k,k1) = Y(k1,1);
        z2(k,k1)  = cos(Y(k1,1))*Ri;
        tturn(k)  = TE; 
        v2(k,k1)  = Ri*Y(k1,2)*sin(Y(k1,1));
        FT2(k,k1)  = m*g+m*(Ri*Y(k1,2)).^2*cos(phi(k,k1))/Ri;
     end
end

 
%% 
% Graphische Ausgabe

% Zugkraft, Geschwindigkeit im Umkehrbereich über t
fig = figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
yyaxis left;
hold on; 
for k=1:4 
    rp(k) = plot(td2(k,:)*1000,FT2(k,:));
    set(rp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
end
grid on
xlabel('Umkehrzeit t in ms','FontSize',14)
ylabel('F_T in N ','FontSize',14)
legend(rp, Para_Str,'location','northwest','numcolumns',2);
legend box off
set(gca,'FontSize',16);

yyaxis right;
hold on; 
for k=1:4 
    rp(k) = plot(td2(k,:)*1000,v2(k,:));
    set(rp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(3));
end
grid on
ylabel('Geschwindigkeit v_z in m/sF_T in N ','FontSize',14)
xlabel('Umkehrzeit t in ms','FontSize',14)
h=title('Zugkraft, Geschwindigkeit im Umkehrbereich');
set(h,'FontSize',14,'FontWeight','normal'); 
legend(rp, Para_Str,'location','northwest','numcolumns',2);
legend box off
set(gca,'FontSize',16);

%%
% Graphische Ausgabe

% Abrollbereich
fig = figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
% Geschwindigkeit über z
subplot(1,2,1);
rp = plot(zspan,v);
for k=1:4 
    set(rp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
end
grid on
xlabel('Fallhöhe h  in m','FontSize',14)
ylabel('v in m/s','FontSize',14)
legend(rp, Para_Str,'location','northwest','numcolumns',2);
legend box off
h=title('Geschwindigkeit');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


% Zugkraft über z
subplot(1,2,2);
hold on
yyaxis left;
lp = plot(zspan,FT);
for k=1:4 
    set(lp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
end
xlabel('Fallhöhe h  in m','FontSize',14)
ylabel('F_T in N','FontSize',14)
set(gca,'FontSize',16);
% Beschleunigung über z
yyaxis right
rp = plot(zspan,a);
for k=1:4 
    set(rp(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(2));
end
xlabel('Fallhöhe h  in m','FontSize',14)
ylabel('a in m/s²','FontSize',14)
grid on
h=title('Zugkraft und a');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);


% Fallhöhe über t
figure()
hold on
for k=1:4
   p(k)=plot(td(k,:), z(k,:));
   set(p(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
   rp1(k) = plot(tfall(k)+td2(k,:),L+z2(k,:));
   set(rp1(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
   rp2(k) = plot(tfall(k)+tturn(k)+td(k,:),zr(k,:));
   set(rp2(k),'Color',Colors(k+1,:), 'LineWidth',2,'LineStyle',Style(1));
end
ylabel('Fallhöhe h  in m','FontSize',14)
xlabel('t in s','FontSize',14)
axis ij
grid on
h=title('z(t)');
legend(p, Para_Str,'location','southwest','numcolumns',1);
legend box off
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




%% Funktionen
% Jo-Jo in Fallphase
function dY = dgl_JoJo1(t,Y,P1)
% Y(1)-  z(t), 
% Y(2)-  Geschwindigkeit dz(t)
rx    = sqrt(P1.Ri^2+P1.d*(P1.L-Y(1))/pi);
a     = 2*P1.g*rx.^2./(P1.Ra^2+2*rx.^2);
dY    = [Y(2);...    
         a];
end


% Jo-Jo in Umkehrphase
function dY = dgl_JoJo2(t,Y,P1)
% Y(1)- Winkel phi(t), 
% Y(2)- Winkelgeschwindigkeit omega(t)
dY    = [Y(2);...    
         -2*sin(Y(1))*P1.g*P1.Ri/P1.Ra/P1.Ra];
end

 
% Ereignisfunktion

function [value,isterminal,direction] = MyEvent1(t,Y,P1)
% Ereignisfunktion bis Eintreten des Gleitens bei Slip Winkel thetaS bzw.
% thetaC
    value = (P1.L-Y(1));    % detect length of line
    isterminal = 1;      % stop the integration
    direction = 0;       % negative direction
end

function [value,isterminal,direction] = MyEvent2(t,Y,P1)
% Ereignisfunktion bis Eintreten des Gleitens bei Slip Winkel thetaS bzw.
% thetaC
    value = (-pi/2-Y(1)); % detect angle pi/2
    isterminal = 1;      % stop the integration
    direction = 0;       % negative direction
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------



