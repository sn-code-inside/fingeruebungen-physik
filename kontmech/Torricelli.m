% -------------------------------------------------------------------------
% Torricelli.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet das Stömungsfeld für den Fluss durch einen Tank
% bei vorausgesetztem Gesetz von Torricelli.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
Nit = 4000;                           % (maximale) Zahl an Iterationen
s = 0.1;                              % Konvergenzparameter

dh = 4.e-3;                           % Abstand zwischen zwei Gitterpunkten
                                      % in m
Nx =  50;                             % Gitterpunkte in x-Richtung
Ny = 200;                             % Gitterpunkte in y-Richtung
NL =  40;                             % Gitterpunkt nx, an dem das
                                      % Ausflussloch beginnt
Nu = 20;                              % Gitterpunkte unter dem Ausfluss

g   = 9.81;                           % Schwerebeschleunigung in m/s^2
nu  = 2.5e-3;                         % kinematische Viskosität

% Setze das Gitter
for i = 1:2*Nx+1
    ValX(i) = dh*(i-1);
end
for i = 1:Ny+1
    ValY(i) = dh*(i-1);
end

% Berechnen der Geschwindigkeiten am unteren und oberen Rand
H = Ny*dh;                            % Höhe der Wassersäule im Tank
vu = sqrt(2*g*H);                     % unterer Rand

%% Bereite die Gitter vor
u = zeros(Nx+1,Ny+1);
Omega = zeros(Nx+1,Ny+1);

%% Integration der Lösung
it = 0;
while it < Nit
    % Interation für u
    Y          = SetzeRaender(u,Omega,g,dh,Nx,Ny,NL,Nu);
    u(:,:)     = Y(1,:,:);
    Omega(:,:) = Y(2,:,:);
    for i = 2:Nx
        for j = 2:Ny
            if (j<=Nu)
                if (i>NL)
                    r1 = s*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) ...
                        +2*dh*dh*Omega(i,j))*0.25-(u(i,j)));
                    u(i,j) = u(i,j) + r1;
                end
            end
            if (j>Nu)
                r1 = s*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) ...
                    +2*dh*dh*Omega(i,j))*0.25-u(i,j));
                u(i,j) = u(i,j) + r1;
            end
        end
    end

    % Iteration für Omega
    Y          = SetzeRaender(u,Omega,g,dh,Nx,Ny,NL,Nu);
    u(:,:)     = Y(1,:,:);
    Omega(:,:) = Y(2,:,:);
    for i = 2:Nx
        for j = 2:Ny
            if (j<=Nu)
                if (i>NL)
                    a1 = Omega(i+1,j)+Omega(i-1,j)+Omega(i,j+1) ...
                        +Omega(i,j-1);
                    a2 = (u(i,j+1)-u(i,j-1))*(Omega(i+1,j)-Omega(i-1,j));
                    a3 = (u(i+1,j)-u(i-1,j))*(Omega(i,j+1)-Omega(i,j-1));
                    r2 = s*(0.25*(a1+(0.25/nu)*(a3-a2))-Omega(i,j));
                    Omega(i,j) = Omega(i,j)+r2;
                end
            end
            if (j>Nu)
                a1 = Omega(i+1,j)+Omega(i-1,j)+Omega(i,j+1)+Omega(i,j-1);
                a2 = (u(i,j+1)-u(i,j-1))*(Omega(i+1,j)-Omega(i-1,j));
                a3 = (u(i+1,j)-u(i-1,j))*(Omega(i,j+1)-Omega(i,j-1));
                r2 = s*(0.25*(a1+(0.25/nu)*(a3-a2))-Omega(i,j));
                Omega(i,j) = Omega(i,j)+r2;
            end
        end
    end

    it = it+1;
end

%% Darstellung der Ergebnisse
% Füge gespiegelte zweite Hälfte entlang der x-Achse für die
% grafische Darstellung hinzu.
u2(1:Nx+1,:)     = u(:,:);
Omega2(1:Nx+1,:) = Omega(:,:);
for i=1:Nx
    u2(i+Nx+1,:)     = u(Nx+1-i,:);
    Omega2(i+Nx+1,:) = Omega(Nx+1-i,:);
end

% Berechne Geschwindigkeitsvektoren
[gY,gX] = gradient(u,dh);

% Selektiere Punkte, die verwendet werden sollen
Dx  = 5;   % jeder Dx-te Punkt entlang der x-Achse wird verwendet
Dy  = 10;  % jeder Dy-te Punkt entlang der y-Achse wird verwendet
NNx = fix(Nx/Dx);
NNy = fix(Ny/Dy);
for i = 1:2*NNx
    ValXt(i) = ValX(Dx*(i-1)+1);
end
for j = 1:NNy         
    ValYt(j) = ValY(Dy*(j-1)+1);
end
for i = 1:NNx
    for j = 1:NNy         
        gXt(i,j) = gX(Dx*(i-1)+1,Dy*(j-1)+1);
        gYt(i,j) = gY(Dx*(i-1)+1,Dy*(j-1)+1);
    end
end
% Erweitere Matrizen um zweite Hälfte des Tanks
gXt2(1:NNx,:) = gXt(1:NNx,:);
gYt2(1:NNx,:) = gYt(1:NNx,:);
for i = 1:NNx
    gXt2(NNx+i,:)   = gXt(NNx-i+1,:);
    gYt2(NNx+i,:)   = -gYt(NNx-i+1,:);
end

% u-Konturplot
figure
subplot(1,3,1)
contour(ValX,ValY,transpose(u2),30,'linewidth',2, 'color', Colors(3,:));
grid on;
xlabel('{\itx}/m');
ylabel('{\ity}/m');
axis([0 0.4 0 0.8])
axis equal
title('{\itu}({\itx},{\ity})','FontSize',14,'FontName','Times', ...
      'FontWeight','normal');
set(gca,'FontSize',16,'FontName','Times');

% Omega-Konturplot (Vortexlinien)
subplot(1,3,2)
contour(ValX,ValY,transpose(Omega2),100,'linewidth',2, 'color', ...
    Colors(3,:));
grid on;
xlabel('{\itx}/m');
ylabel('{\ity}/m');
axis([0 0.4 0 0.8])
axis equal
title('{\Omega}({\itx},{\ity})','FontSize',16,'FontName','Times', ...
      'FontWeight','normal');
set(gca,'FontSize',16,'FontName','Times');

subplot(1,3,3)
axis([0 0.4 0 0.8])
axis equal
hold on
q= quiver(ValXt,ValYt,transpose(gYt2),-transpose(gXt2),'linewidth',2, ...
    'color', Colors(3,:));
grid on;
xlabel('{\itx}/m');
ylabel('{\ity}/m');
axis equal
axis([0 0.4 0 0.8])
title('{\it\bfv}({\itx},{\ity})','FontSize',16,'FontName','Times', ...
      'FontWeight','normal');
set(gca,'FontSize',16,'FontName','Times');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Funktion
% Ränder setzen
function Y = SetzeRaender(u,Omega,g,dh,Nx,Ny,NL,Nu)
  % unter dem Loch  
  for i = NL:Nx+1
      u(i,1)       = u(i-1,2);       % du/dy = vx = 0
      Omega(i-1,1) = Omega(i-1,2);   % Wasser am Boden
      for j = 1:Nu+1
          if (i==NL)
              vy = 0;
          end
          if (i==Nx+1)
              vy = -sqrt(2.*g*dh*(Ny+Nu-j));
          end
          if (i==Nx)
              vy = -0.5*sqrt(2.*g*dh*(Ny+Nu-j));
          end
          u(i,j) = u(i-1,j)-vy*dh;    % du/dx = -vy
      end
  end
  
  % rechter Rand
  for j = 2:Ny+1
      vy          = sqrt(2.*g*dh*(Ny+1-j));
      u(Nx+1,j)   = u(Nx,j)-vy*dh;       % du/dx = -vy
      u(Nx+1,j)   = u(Nx+1,j-1);         % du/dy = vx = 0
      Omega(Nx,j) = -2*(u(Nx,j)-u(Nx,j-1))/dh^2;
  end
  
  % Boden, links vom Loch
  for i = 1:NL-1
      u(i,Nu)     = u(i,Nu-1);          % du/dy = vx = 0
      Omega(i,Nu) = -2*(u(i,1)-u(i,2));
  end
  
  % oberer Rand
  for i = 1:Nx
      u(i,Ny+1)   = u(i,Ny);            % du/dy = vx = 0
      Omega(i,Ny) = 0;
  end
    
  % linker Rand
  for j = Nu:Ny
      Omega(1,j) = -2*(u(1,j)-u(2,j))/dh^2;
      u(1,j)     = u(2,j);              % du/dx = -vy = 0
  end
  
  Y(1,:,:) = u(:,:);
  Y(2,:,:) = Omega(:,:);
end