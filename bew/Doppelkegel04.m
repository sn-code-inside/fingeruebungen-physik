% -------------------------------------------------------------------------
% Doppelkegel04.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Aufwärtsrollender Doppelkegel
% 
% Programm berechnet Lösungen den Phasenraum aus den
% Lagrange-Gleichungen des aufwärtsrollenden Doppelkegels und die
% Energiekonversion für verschiedene Parameter
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
H       = 6.5;              % Länge Schaukel cm
R       = 3.0;              % Länge Oberkörper cm 
phi     = deg2rad(15.3);    % Öffnungswinkel Schienen
theta   = deg2rad(06.5);    % Anstiegswinkel Schienen
psi     = atan(R/H);        % Öffnungswinkel Kegel
M       = 122.5;            % Masse Doppelkegel in g
J       = 3*M*R^2/10;       % Trägheitsmoment Doppelkegel in g cm²
g       = 981;              % g in cm/s

alpha   = asin(tan(psi)*tan(phi));  %Richtungswinkel alpha
q0      = 0;                % Anfangsgposition gen. KO des Doppelkegels
dq0     = 0;                % Anfangsgeschwindigkeit gen. KO des Doppelkegels
ys0     = H*tan(psi)*sin(alpha);
                             % Anfangsgposition ys im KOS Sigma


fprintf('\n ');
fprintf('\n phi   = %4.2f°', rad2deg(phi));
fprintf('\n psi   = %4.2f°', rad2deg(psi));
fprintf('\n theta = % 4.2f°', rad2deg(theta));
fprintf('\n alpha = % 4.2f°', rad2deg(alpha));
fprintf('\n ');
fprintf('\n H  = %8.2f cm', H);
fprintf('\n R  = %8.2f cm', R);
fprintf('\n M  = %8.2f g', M);
fprintf('\n J  = %8.2f gcm²', J);
fprintf('\n ');
fprintf('\n ys0= %8.2f cm', ys0);
fprintf('\n g  = %8.2f cm/s', g);
fprintf('\n ');
fprintf('\n Rollbedingung erfüllt ? ');
b1 = tan(theta);
b2 = tan(phi)*tan(psi)/sqrt(1-(tan(phi)*tan(psi))^2);
if b1 < b2
    fprintf('\n %s  < %s ! Rollbedingung erfüllt! ',...
        num2str(b1,4), num2str(b2,4));
    fprintf('\n ');
else
    fprintf('\n %s  > %s ! Rollbedingung nicht (!) erfüllt! ',...
        num2str(b1,4), num2str(b2,4));
    fprintf('\n ');
end

fprintf('\n ');

%% Berechnungen

% Anfangswerte
iplotE =3;

for k=1:iplotE  
    phiv     =[deg2rad(15.5),deg2rad(15.5),deg2rad(15.0),deg2rad(16.0)] ; 
    % Öffnungswinkel Schienen 
    thetav   =[deg2rad(06.5),deg2rad(5.5),deg2rad(5.5),deg2rad(06.5)];  
    % Anstiegswinkel Schienen 
    Jv       =[330.75,330.75,165.75,165.75];           
    % Trägheitsmoment in g cm²
    phi     = phiv(k);
    theta   = thetav(k);
    alpha   = asin(tan(psi)*tan(phi));  % Hilfswinkel
    J       = Jv(k);
    qmax(k)  = 0.9999*H*tan(psi)/tan(alpha);

    q(k,:) = linspace(0,qmax(k),100);
    E0      = 0;
    dq(k,:) = sqrt(2*(E0+M*g*q(k,:)*sin(alpha-theta))./...
                  (M+(J./(H*tan(psi)-q(k,:)*tan(alpha)).^2)));
    C3   = H*tan(psi)/tan(alpha);

    U(k,:)      = E0+M*g*q(k,:)*sin(alpha-theta)/1e5;
    T_trans(k,:)= M*(tan(alpha))^2.*(C3-q(k,:)).^2.*(E0+M*g*q(k,:)*sin(alpha-theta))./...
                 (M*(tan(alpha))^2.*(C3-q(k,:)).^2 +J)/1e5;
    T_rot(k,:)  = J*(E0+M*g*q(k,:)*sin(alpha-theta))./...
                   (M*(tan(alpha))^2.*(C3-q(k,:)).^2 +J)/1e5;
end
%% Graphische Ausgabe
% Lösung Phasenraum
figure();
hold on
for k=1:iplotE
    lp(k)=plot(q(k,:),dq(k,:),'Color',Colors(k,:), 'LineWidth',2);
    line([qmax(k) qmax(k)],[0, 1.0*max(dq(k,:))],'Color',Colors(8,:),...
         'LineWidth',1,'LineStyle',Style(3));
end
axis([0 25 0 1.0*max(max(dq(:,:)))]);
grid on
xlabel('q in cm','FontSize',14)
ylabel('dq in cm/s','FontSize',14)
legend(lp,'Parameter A','Parameter B','Parameter C','location','south');
legend box off
h=title('Phasenraum');
set(h,'FontSize',12,'FontWeight','normal'); 
set(gca,'FontSize',16);

% Energiegleichung
figure();
hold on
rp(1)=plot(q(1,:),T_trans(1,:));
rp(2)=plot(q(1,:),T_rot(1,:));
rp(3)=plot(q(1,:),U(1,:));
for k=1:3 
    set(rp(k),'Color',Colors(k+1,:),'LineWidth',2,'LineStyle',Style(1));
    line([qmax(1) qmax(1)],[0, max(E0+U(1,:))],'Color',Colors(4,:),'LineWidth',1);
end
text(20,1.1*max(U(1,:)),'Parameter A','Color',Colors(4,:));
rp(4)=plot(q(3,:),T_trans(3,:));
rp(5)=plot(q(3,:),T_rot(3,:));
rp(6)=plot(q(3,:),U(3,:));
for k=1:3 
    set(rp(k+3),'Color',Colors(k+1,:),'LineWidth',2,'LineStyle',Style(3));
    line([qmax(3) qmax(3)],[0, max(E0+U(3,:))],'Color',Colors(3,:),...
         'LineWidth',1,'LineStyle',Style(3))
end
text(20,1.1*max(U(3,:)),'Parameter C','Color',Colors(4,:));
axis([0 25 0 round(1.2*max(E0+U(3,:))*10)/10]);
grid on
xlabel('q in cm','FontSize',14)
ylabel('E(q) in Ws ','FontSize',14)
legend('T_{trans}', 'T_{rot}', 'U', 'location','northwest');
legend box off
h=title('Energieverteilung');
set(h,'FontSize',12,'FontWeight','normal');
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


