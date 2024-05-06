% -------------------------------------------------------------------------
% Simuslauf.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Sinuslauf bei starren Radachsen auf Schienen
% 
% Programm berechnet Frequenz und Stabilität des 
% Sinuslaufs bei starren Radachsen auf Schienen
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];


%% Numerik
% Parameter

H       = 1.435/2;         % Schienenbreite in m
R0      = 0.4875;          % Rollradius in m
psi     = deg2rad(1/40);   % Konizität
v0      = 50.0;            % Geschwindigkeit in m/s
g       = 9.81;            % in m/s^2
C       = 1;               % dimensionslose Konstante (C= Jy/mH^2)

fprintf('\n ');
fprintf('\n Klingel-Frequenz Parameter:');
fprintf('\n ');
fprintf('\n Schienenbreite :             2H   = %8.3f m', 2*H);
fprintf('\n Rollradius :                 R0   = %8.3f m', R0);
fprintf('\n Konizität :                  psi  = %8.3f ° ', rad2deg(psi));
fprintf('\n Geschwindigkeit :            v0   = %8.3f m/s', v0);
fprintf('\n ');


%% Berechnungen Einfache Klingel-Formel

psiv = linspace(0.1,2,100);
psiv = deg2rad(psiv);
Hv     = [H/2,H,2*H];
Cv     = [1,1,1];

for k= 1:length(Hv)
  fK(k,:) = v0.*sqrt(psiv(:)/Hv(k)/R0)/2/pi;
  lambdaK(k,:) = 2*pi*sqrt(R0*Hv(k)./psiv);
  parastr(k,:)  = string(strcat(strcat('H = ',num2str(Hv(k)*1000,4)),' mm'));
end
% Graphische Ausgabe

% Frequenz 
figure();
for k= 1:length(Hv)
  p1(k) = plot(rad2deg(psiv(:)), fK(k,:),'Color',Colors(k,:),...
         'LineWidth',2);
  hold on
end
grid on
xlabel('\psi in °','FontSize',14)
ylabel('Frequenz f_K in 1/s','FontSize',14)
h=title('Klingel-Formel ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p1,parastr(:,:), 'location','northwest',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);

% Wellenlänge 
figure();
for k= 1:length(Hv)
  p2(k) = plot(rad2deg(psiv(:)), lambdaK(k,:),'Color',Colors(k,:),...
         'LineWidth',2,'LineStyle',Style(1));
  hold on;
end
grid on
xlabel('\psi in °','FontSize',14)
ylabel('\lambda_K in m','FontSize',14)
h=title('Klingel-Formel');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p2,parastr(:,:), 'location','northeast',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);


%% Dynamische Berechnungen Sinuslauf-Frequenz

m      = 1000;           % in kg
Jyv    = 2.5*m*Hv.^2;    % in kgm^2
omega0 = v0/R0;

for k= 1:length(Hv)
  fS(k,:) = sqrt(2*m*g./(Jyv(k)*psiv(:)/Hv(k)-m*R0+m*Hv(k)./psiv(:)))/2/pi;
end
% Graphische Ausgabe

% Frequenz 
figure();
for k= 1:length(Hv)
  p1(k) = plot(rad2deg(psiv(:)), fS(k,:),'Color',Colors(k,:),...
         'LineWidth',2);
  hold on
end
grid on
xlabel('\psi in °','FontSize',14)
ylabel('Frequenz f_S in 1/s','FontSize',14)
h=title('Sinuslauf ');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p1,parastr(:,:), 'location','northwest',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);


%% Eigenwerte und Stabilitaetsdiagramm

eta    = 4e3;            % in kg/s
etav   = [eta/2, eta, eta*1.5];    
m      = 25000;          % in kg
Jz     = 2.27e6;         % in kgm^2
Jy     = 0.26e6;         % in kgm^2
omega0 = v0/R0;

fprintf('\n ');
fprintf('\n Sinuslauf- Dynamische Betrachtung Parameter:');
fprintf('\n ');
fprintf('\n Schienenbreite :             2H   = %8.3f m', 2*H);
fprintf('\n Rollradius :                 R0   = %8.3f m', R0);
fprintf('\n Konizität :                  psi  = %8.3f ° ', rad2deg(psi));
fprintf('\n Geschwindigkeit :            v0   = %8.3f m/s', v0);
fprintf('\n Masse :                      m    = %8.3f kg', m);
fprintf('\n Trägheitsmoment :            Jz   = %8.3e kg m²', Jz);
fprintf('\n Trägheitsmoment :            Jy   = %8.3e kg m²', Jy);
fprintf('\n Reibungskoeefizient :        eta  = %8.3e kg/s', eta);
fprintf('\n ');

psiv = linspace(0.01,0.06,100);

for k2 = 1:length(etav)
    eta = etav(k2);
    parastr2(k2,:)  = string(strcat(strcat('\eta = ',...
                       num2str(etav(k2),4)),' kg/s'));
    Q = - Jy*psiv/H - m*R0 + m*H./psiv;
    M11 = 2*m*g./Q;
    M12 = -2*eta*omega0*R0./psiv./Q;
    M21 = 2*eta*omega0*H^2*psiv/Jz;
    M22 = -2*m*g*H*psiv/Jz;
    for k=1:length(psiv)   
      M = [M11(k),M12(k);M21(k),M22(k)];
      lambda(k,:) = eig(M);
      for k1 = 1: length(lambda(k,:))
          if imag(lambda(k,k1)) ~= 0 
             lambda(k,k1) = NaN;
          end
      end
    end
    lambda_x(:,k2) = lambda(:,1);
    lambda_Y(:,k2) = lambda(:,2);
    Trace  = (M11+ M22);
    Det    = (M11.*M22-M12.*M21);
    LAMBDA(k2,:) = Det(:)./Trace(:).^2;
end



% Eigenwerte 
figure();
for k2=1:length(etav)
    p3(k2) = plot(rad2deg(psiv), LAMBDA(k2,:),'Color',Colors(k2+1,:),...
             'LineWidth',2);
    hold on
    p4(k2) = plot(rad2deg(psiv), lambda_x(:,k2),'Color',Colors(k2+1,:),...
         'LineWidth',1,'LineStyle',Style(2));
    p5(k2) = plot(rad2deg(psiv), lambda_Y(:,k2),'Color',Colors(k2+1,:),...
         'LineWidth',1,'LineStyle',Style(3));
% 
end
p3(k2+1) = line([0,rad2deg(psiv(end))],[0.25 0.25],'Color',Colors(8,:),...
           'LineWidth',2);
parastr2(k2+1,:) = ' Stabilitätsgrenze \Lambda';
grid on
xlabel('\psi in °','FontSize',14)
ylabel('\Lambda, \lambda in 1/s²','FontSize',14)
axis([0 rad2deg(psiv(end)) -0.1 0.8]);
h=title('Sinuslauf-Stabilität');
set(h,'FontSize',12,'FontWeight','normal'); 
legend(p3,parastr2(:,:), 'location','northeast',...
         'NumColumns',1);
legend box off
set(gca,'FontSize',16);
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


