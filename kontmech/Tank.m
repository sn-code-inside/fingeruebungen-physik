% -------------------------------------------------------------------------
% Tank.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik des Kontinuums" aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Programm berechnet das Stömungsfeld für den Fluss durch einen voll
% gefüllten Tank, bei dem links und rechts zwei gleich große Löcher
% einen Durchfluss von links nach rechts ermöglichen. Die
% Anfangsbedingungen sind so gewählt, dass ein Fluss am oberen Rand
% des Tanks beginnt.
% -------------------------------------------------------------------------

clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

%% Parameter
Nit = 100000;                       % (maximale) Zahl an Iterationen
s = 0.1;                            % Konvergenzparameter

dh = 3.e-3;                         % Abstand zwischen zwei Gitterpunkten
                                    % in m
Nx = 100;                           % Gitterpunkte in x-Richtung
Ny = 100;                           % Gitterpunkte in y-Richtung

Nlu = 60;                           % unterer Rand des Einflusslochs
Nlo = 70;                           % oberer Rand des Einflusslochs
Nru = 20;                           % unterer Rand des Ausflusslochs
Nro = 30;                           % oberer Rand des Ausflusslochs
                                      
nu  = 2.5e-3;                       % kinematische Viskosität

% Setze das Gitter
for i = 1:Nx+1
    ValX(i) = dh*(i-1);
end
for i = 1:Ny+1
    ValY(i) = dh*(i-1);
end

% Geschwindigkeiten
v0 = 0.18;                             % Fluss durch die Löcher
vo = 0.3;                              % Geschwindigkeit am oberen Rand

%% Bereite die Gitter vor
u = zeros(Nx+1,Ny+1);
Omega = zeros(Nx+1,Ny+1);


fprintf('\n');
fprintf('\n Numerische Integration läuft');
fprintf('\n Bitte etwas Geduld ...');
fprintf('\n');

%% Integration der Lösung
it = 0;
while it < Nit
    % Interation für u
    Y          = SetzeRaender(u,Omega,dh,Nx,Ny,Nlu,Nlo,Nru,Nro,v0,vo);
    u(:,:)     = Y(1,:,:);
    Omega(:,:) = Y(2,:,:);
    for i = 2:Nx
        for j = 2:Ny
            r1 = s*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1) ...
                +2*dh*dh*Omega(i,j))*0.25-u(i,j));
            u(i,j) = u(i,j) + r1;
        end
    end

    % Iteration für Omega
    Y          = SetzeRaender(u,Omega,dh,Nx,Ny,Nlu,Nlo,Nru,Nro,v0,vo);
    u(:,:)     = Y(1,:,:);
    Omega(:,:) = Y(2,:,:);
    for i = 2:Nx
        for j = 2:Ny
            a1 = Omega(i+1,j)+Omega(i-1,j)+Omega(i,j+1)+Omega(i,j-1);
            a2 = (u(i,j+1)-u(i,j-1))*(Omega(i+1,j)-Omega(i-1,j));
            a3 = (u(i+1,j)-u(i-1,j))*(Omega(i,j+1)-Omega(i,j-1));
            r2 = s*(0.25*(a1+(0.25/nu)*(a3-a2))-Omega(i,j));
            Omega(i,j) = Omega(i,j)+r2;
        end
    end

    it = it+1;
end

fprintf('\n');
fprintf('\n Numerische Integration fertig');
fprintf('\n');
pause(0.2);


%% Darstellung der Ergebnisse
% Berechne Geschwindigkeitsvektorenm
[gY,gX] = gradient(u,dh);

% Selektiere Punkte, die verwendet werden sollen
Dx  = 5;  % jeder Dx-te Punkt entlang der x-Achse wird verwendet
Dy  = 5;  % jeder Dy-te Punkt entlang der y-Achse wird verwendet
NNx = fix(Nx/Dx);
NNy = fix(Ny/Dy);
for i = 1:NNx
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


% u-Konturplot (Flusslinien)
figure
subplot(1,3,1)
contour(ValX,ValY,transpose(u),30,'linewidth',2, 'color', Colors(3,:));
grid on;
xlabel('{\itx}/m');
ylabel('{\ity}/m');
axis([0 0.3 0 0.3])
axis equal
title('{\itu}({\itx},{\ity})','FontSize',14,'FontName','Times', ...
      'FontWeight','normal');
set(gca,'FontSize',16,'FontName','Times');

% Omega-Konturplot (Vortexlinien)
subplot(1,3,2)
contour(ValX,ValY,transpose(Omega),100,'linewidth',2, 'color', Colors(3,:));
grid on;
xlabel('{\itx}/m');
ylabel('{\ity}/m');
axis([0 0.3 0 0.3])
axis equal
title('{\Omega}({\itx},{\ity})','FontSize',14,'FontName','Times', ...
      'FontWeight','normal');
set(gca,'FontSize',16,'FontName','Times');

subplot(1,3,3)
axis([0 0.3 0 0.3])
axis square
hold on
q= quiver(ValXt,ValYt,transpose(gYt),-transpose(gXt),'linewidth',2, ...
    'color', Colors(3,:));
grid on;
xlabel('{\itx}/m');
ylabel('{\ity}/m');
axis([0 0.3 0 0.3])
axis square
title('{\bfv}({\itx},{\ity})','FontSize',14,'FontName','Times', ...
      'FontWeight','normal');
set(gca,'FontSize',16,'FontName','Times');

% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------



%% Funktion
% Ränder setzen
function Y = SetzeRaender(u,Omega,dh,Nx,Ny,Nlu,Nlo,Nru,Nro,v0,vo)
  % oberer Rand
  for i = 2:Nx+1
      u(i-1,Ny+1)   = u(i,Ny+1);                  % du/dx = -vy = 0
      Omega(i,Ny+1) = 2*(u(i,Ny+1)-u(i,Ny))/dh^2 - 2*vo/dh;
  end

  % unterer Rand
  for i = 2:Nx+1
      Omega(i,1) = 2*(u(i,1)-u(i,2))/dh^2;
      u(i,2)     = u(i,1);                  % du/dy = vx = 0
      u(i-1,1)   = u(i,1);                  % du/dx = -vy = 0
  end

 % linker Rand
 for j = 2:Ny+1
     if (j<Nlu)
         % unter dem Einflussloch
         Omega(1,j) = 2*(u(1,j)-u(2,j))/dh^2;
         u(2,j)     = u(1,j);                 % du/dx = -vy = 0
         u(1,j)     = u(1,j-1);               % du/dy = vx = 0
     elseif ((j>=Nlu) && (j<=Nlo))
         % entlang des Einflusslochs
         Omega(2,j) = Omega(1,j);  
         u(1,j)     = u(1,j-1) + v0*dh;       % du/dy = vx = v0
     elseif (j>Nlo)
         % oberhalb des Einflusslochs
         u(1,j)     = u(1,j-1);               % du/dy = vx = 0
         u(2,j)     = u(1,j);                 % du/dx = -vy = 0 
         Omega(1,j) = 2*(u(1,j)-u(2,j))/dh^2;
     end
 end

 % rechter Rand
 for j = 2:Ny+1
     if (j<Nru)
         % unter dem Ausflussloch
         Omega(Nx,j) = 2*(u(Nx,j)-u(Nx+1,j))/dh^2;
         u(Nx,j)     = u(Nx+1,j);             % du/dx = -vy = 0
         u(Nx+1,j)   = u(Nx+1,j-1);           % du/dy = vx = 0
     elseif ((j>=Nru) && (j<=Nro))
         % entlang des Ausflusslochs
         u(Nx+1,j)   = u(Nx,j);               % du/dx = -vy = 0
         u(Nx+1,j)   = u(Nx+1,j-1) + v0*dh;   % du/dy = vx = v0
         Omega(Nx,j) = Omega(Nx+1,j);
     elseif (j>Nro)
         % oberhalb des Ausflusslochs
         u(Nx+1,j)   = u(Nx+1,j-1);           % du/dy = vx = 0
         u(Nx,j)     = u(Nx+1,j);             % du/dx = -vy = 0 
         Omega(Nx,j) = 2*(u(Nx,j)-u(Nx+1,j))/dh^2;
     end
 end
  
  Y(1,:,:) = u(:,:);
  Y(2,:,:) = Omega(:,:);
end