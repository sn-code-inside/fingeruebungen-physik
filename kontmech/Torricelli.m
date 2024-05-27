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
Nit = 501;                            % (maximale) Zahl an Iterationen
s = 0.1;                              % Konvergenzparameter

dh = 4.e-3;                           % Abstand zwischen zwei Gitterpunkten
                                      % in m
Nx = 200;                             % Gitterpunkte in x-Richtung
Ny = 200;                             % Gitterpunkte in y-Richtung
NL = 180;                             % Gitterpunkt nx, an dem das
                                      % Ausflussloch beginnt
Nu = 20;                              % Gitterpunkte unter dem Ausfluss

g = 9.81;                             % Schwerebeschleunigung in m/s^2
nu = 0.5;                             % kinematic viscosity

% Berechnen der Geschwindigkeiten am unteren und oberen Rand
H = Ny*dh;                            % Höhe der Wassersäule im Tank
vu = sqrt(2*g*H);                     % unterer Rand
vo = vu*(Ny-NL)/Ny;                   % oberer Rand aus Massenerhaltung

% Berechnen der Reynoldszahl
R = vo*dh/nu;

%% Bereite die Gitter vor
u = zeros(Nx+1,Ny+1);
Omega = zeros(Nx+1,Ny+1);

%% Integration der Lösung
it = 0;
while it < Nit
    % Interation für u
    SetzeRaender(u,Omega,g,dh,Nx,Ny,NL,Nu);
    for i = 2:Nx
        for j = 2:Ny
            if (j<=Nu)
                if (i>NL)
                    r1 = s*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) ...
                        +dh*dh*Omega(i,j))*0.25-(u(i,j)));
                    u(i,j) = u(i,j) + r1;
                end
            end
            if (j>Nu)
                r1 = s*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) ...
                    +dh*dh*Omega(i,j))*0.25-u(i,j));
                u(i,j) = u(i,j) + r1;
            end
        end
    end

    SetzeRaender(u,Omega,g,dh,Nx,Ny,NL,Nu);
    % Iteration für w
    for i = 2:Nx
        for j = 2:Ny
            if (j<=Nu)
                if (i>NL)
                    a1 = Omega(i+1,j)+Omega(i-1,j)+Omega(i,j+1) ...
                        +Omega(i,j-1);
                    a2 = (u(i,j+1)-u(i,j-1))*(Omega(i+1,j)-Omega(i-1,j));
                    a3 = (u(i+1,j)-u(i-1,j))*(Omega(i,j+1)-Omega(i,j-1));
                    r2 = s*((0.25*(a1+0.25*R)*(a3-a2))-Omega(i,j));
                    Omega(i,j) = Omega(i,j)+r2;
                end
            end
            if (j>Nu)
                a1 = Omega(i+1,j)+Omega(i-1,j)+Omega(i,j+1)+Omega(i,j-1);
                a2 = (u(i,j+1)-u(i,j-1))*(Omega(i+1,j)-Omega(i-1,j));
                a3 = (u(i+1,j)-u(i-1,j))*(Omega(i,j+1)-Omega(i,j-1));
                r2 = s*((0.25*(a1+0.25*R)*(a3-a2))-Omega(i,j));
                Omega(i,j) = Omega(i,j)+r2;
            end
        end
    end

    it = it+1;
end

%% Darstellung der Ergebnisse

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------


%% Funktion
% Raender setzen
function SetzeRaender(u,Omega,g,dh,Nx,Ny,NL,Nu)
  % unter dem Loch  
  for i = NL:Nx+1
      u(i,1) = u(i,2);               % du/dy = vx = 0
      Omega(i-1,1) = Omega(i-1,2);   % Wasser am Boden
      for j = 1:Nu
          if (i==NL)
              vy = 0;
          end
          if (i==Nx+1)
              vy = -sqrt(2.*g*dh*(Ny+Nu-j));
          end
          if (i==Nx)
              vy = -0.5*sqrt(2.*g*dh*(Ny+Nu-j));
          end
          u(i,j) = u(i-1,j)-vy*dh;   % du/dx = -vy
      end
  end

  % rechter Rand
  for j = 2:Ny+1
      vy = sqrt(2.*g*dh*(Ny+1-j));
      u(Nx,j) = u(Nx-1,j)-vy*dh;     % du/dx = -vy
      u(Nx,j) = u(Nx,j-1);           % du/dy = vx = 0
      w(Nx,j) = -2*(u(Nx,j)-u(Nx,j-1))/dh^2;
  end

  % Boden, links vom Loch
  for i = 1:NL-1
      u(i,Nu) = u(i,Nu-1);          % du/dy = vx = 0
      w(i,Nu) = -2*(u(i,1)-u(i,2));
  end

  % oberer Rand
  for i = 1:Nx+1
      u(i,Ny) = u(i,Ny-1);          % du/dy = vx = 0
      w(i,Ny) = 0;
  end
  
  % linker Rand
  for j = Nu:Ny+1
      w(1,j) = -2*(u(1,j)-u(2,j))/dh^2;
      u(1,j) = u(2,j);              % du/dx = vy = 0 
  end
end