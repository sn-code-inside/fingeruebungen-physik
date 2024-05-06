% -------------------------------------------------------------------------
% .m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Bungee Jumping Kette
% 
% Programm berechnet Dynamik des Bungee Jumping auf Basis der
% Lagrange-Gleichungen  
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
L       = 20.0;             % Länge Seil m
mS      = 50.0;             % Masse des Seils kg
mK      = 40.0;             % Masse des Körpers des Springers kg
mu      = mS/L;             % Massendichte in kg/m 
gamma   = mS/mK;

g       = 9.81;             % g in m/s^2 
z0      = 0;                % Anfangsgposition der Kette
dz0     = 0;                % Anfangsgeschwindigkeit der Kette
tmax    = sqrt(L/g);        % Fallzeit freier Fall

fprintf('\n ');
fprintf('\n L   = %8.2f m', L);
fprintf('\n mu  = %8.2f kg/m', mu);
fprintf('\n mS  = %8.2f kg/m', mS);
fprintf('\n mK  = %8.2f kg/m', mK);
fprintf('\n ');
fprintf('\n z0  = %8.2f m', z0);
fprintf('\n dz0 = %8.2f m/s', dz0);
fprintf('\n g   = %8.2f m/s^2', g);
fprintf('\n tmax= %8.2f m/s', tmax);
fprintf('\n ');


%% Berechnungen

% Anfangswerte
AB=[z0;dz0]; % AB für ode45


for m= 1:3 
    mK = m*0.5*mS;
    gamma(m) = mK/mS;
    % Numerische Lösung volle Lagrangegleichung
    opt=odeset('AbsTol',1.e-7,'RelTol',1.e-6);
    % Numerische Lösung LGL
    [t,Y]=ode45(@dgl_BungeeJump,[0.0,tmax],AB,opt,g,mK,mS,L); 
    for k=1:length(t)
        if Y(k,1) > L/2
           Y(k,1) = NaN;
           Y(k,2) = NaN;
        end
    end
    
    switch m
    case 1
        t1   = t;
        z1   = Y(:,1);
        dz1  = Y(:,2);
        ddz1 = g*(1+(mS.*z1.*(4*mK*L+2*mS*L-mS.*z1))./(2*(mS*L-mS.*z1+2*mK*L).^2));
   case 2
        t2   = t;
        z2   = Y(:,1);
        dz2  = Y(:,2);
        ddz2 = g*(1+(mS.*z2.*(4*mK*L+2*mS*L-mS.*z2))./(2*(mS*L-mS.*z2+2*mK*L).^2));
    case 3
        t3   = t;
        z3   = Y(:,1);
        dz3  = Y(:,2);
        ddz3 = g*(1+(mS.*z3.*(4*mK*L+2*mS*L-mS.*z3))./(2*(mS*L-mS.*z3+2*mK*L).^2));
    end
end

zf = -0.5*g*t2.^2;
dzf= -g*t2;

%% Graphische Ausgabe

% Lösung Lagrangegleichung z(x)
fig=figure();
set(fig,'defaultAxesColorOrder',[Colors(3,:); Colors(2,:)]);
subplot(1,3,1)

hold on
yyaxis right
plot(t1,abs(z1)-abs(zf),'Color',Colors(2,:), 'LineWidth',2);
plot(t2,abs(z2)-abs(zf),'Color',Colors(2,:), 'LineWidth',2);
% plot(t3,abs(z3)-abs(zf),'Color',Colors(2,:), 'LineWidth',2);
ylabel('Abweichung in m','FontSize',14)
yyaxis left
plot(t2,abs(z2),'Color',Colors(3,:), 'LineWidth',2);
plot(t2,abs(zf),'Color',Colors(4,:), 'LineWidth',2,'LineStyle',Style(3));
axis([0 tmax 0 1.2*L/2]);
grid on
ylabel('z in m','FontSize',14)
xlabel('t in s','FontSize',14)
h=title('Numerische Lösung Lagrange-Gl.');
set(h,'FontSize',14,'FontWeight','normal'); 
set(gca,'FontSize',16);

subplot(1,3,2)
hold on
yyaxis right
plot(z1,abs(dz1)-abs(dzf),'Color',Colors(2,:), 'LineWidth',2);
plot(z2,abs(dz2)-abs(dzf),'Color',Colors(2,:), 'LineWidth',2);
plot(z3,abs(dz3)-abs(dzf),'Color',Colors(2,:), 'LineWidth',2);
ylabel('Abweichung in m/s','FontSize',14)
yyaxis left
plot(z2,abs(dz2),'Color',Colors(3,:), 'LineWidth',2);
plot(z2,abs(dzf),'Color',Colors(4,:),'LineWidth',2,'LineStyle', Style(3));
axis([0 L/2 0 1.2*max(abs(dz2))]);
grid on
ylabel('Geschwindigkeit in m/s','FontSize',14)
xlabel('z in m','FontSize',14)
set(gca,'FontSize',16);

subplot(1,3,3)
hold on
yyaxis right
plot(z1,abs(ddz1)-z1*g./z1,'Color',Colors(2,:), 'LineWidth',2);
plot(z2,abs(ddz2)-z2*g./z2,'Color',Colors(2,:), 'LineWidth',2);
plot(z3,abs(ddz3)-z3*g./z3,'Color',Colors(2,:), 'LineWidth',2);
ylabel('Abweichung in m/s^2','FontSize',14)
yyaxis left
plot(abs(z2),abs(ddz2),'Color',Colors(3,:), 'LineWidth',2);
plot(abs(zf),z2*g./z2,'Color',Colors(4,:),'LineWidth',2,'LineStyle', Style(3));
axis([0 L/2 0 1.2*max(abs(ddz2))]);
grid on
ylabel('Beschleunigung in m/s^2','FontSize',14)
xlabel('z in m','FontSize',14)
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------




%% Funktionen

% Lagrange-Gleichung
function dY = dgl_BungeeJump(t,Y,g,mK,mS,L)
    % Y(1)- Position z(t), 
    % Y(2)- Geschwindigkeit dz(t)
    z = Y(1);
    dY    = [Y(2);...    
         g*(1+(mS*z*(4*mK*L+2*mS*L-mS*z))/(2*(mS*L-mS*z+2*mK*L).^2))];
end
% -------------------------------------------------------------------------
% Ende Funktionen
% -------------------------------------------------------------------------


