% -------------------------------------------------------------------------
% LogistischeGleichung02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel "Physik der Bewegung" aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Getriebener Nichtlinearer Oszillator/Pendel
% Logistische Gleichung - Bifurkationsgraph und Ljapunow-Exponent
% -------------------------------------------------------------------------
clc
clear all
close all 
addpath('./IncludeFolder/')
addpath('./Data/')
Colors = GetColorLines;
Style = ["-", "-.", ":", "--", ":"];

% eine Trajektorie
Gamma0 = 2.5;          
Gamma1 = 4.0;          
NGamma = 2500;      % Varianten von Gamma
Gamma  = Gamma0:(Gamma1-Gamma0)/(NGamma-1):Gamma1;
NPts   = 500;       % Zahl der Zyklen
NPer   = 10;        % Zahl der Bifurkationsperioden 
fder   = @(x) 1-2*x;
for m=1:NGamma
    x(1) = 0.2;
    x(2) = Gamma(m)*x(1)*(1-x(1));
    for k = 2:NPts+1
          x(k+1) = Gamma(m)*x(k)*(1-x(k)); 
    end
    for k = 1:NPer
          xp(k,m) = x(NPts+2-k); 
    end
    lambda(m) = sum(log(abs(Gamma(m)*fder(x))))/(NPts+1);
end

figure();
% Bifurkationsdiagramm
hold on
for k = 1:NPer
      plot(10,10,'b.','MarkerSize',3);
      %plot(Gamma,xp(k,:),'b.','MarkerSize',3);

end
hold off
axis tight
axis([Gamma0 Gamma1 0 1]);
grid on
xlabel('$\Gamma$','interpreter','latex','FontSize',14)
ylabel('$x_{\infty}$','interpreter','latex','FontSize',14)

figure();
% Bifurkationsdiagramm
hold on
for k = 1:NPer
      plot(Gamma,xp(k,:),'b.','MarkerSize',3);

end
hold off
axis tight
axis([Gamma0 Gamma1 0 1]);


figure();
% Bifurkationsdiagramm
hold on
for k = 1:NPer
      plot(Gamma,xp(k,:),'b.','MarkerSize',3);
end
hold off
axis tight
axis([Gamma0 Gamma1 0 1]);
grid on
xlabel('$\Gamma$','interpreter','latex','FontSize',14)
ylabel('$x_{\infty}$','interpreter','latex','FontSize',14)

figure();
% Ljapunow-Exponent
hold on
plot(Gamma,lambda,'color',Colors(2,:),'LineWidth',2);
line([Gamma0 Gamma1], [0 0], 'Color', Colors(3,:),'LineStyle', Style(3),...
       'LineWidth',2);
grid on
axis([Gamma0 Gamma1 -3 1]);
line 
xlabel('$\Gamma$','interpreter','latex','FontSize',14)
ylabel('$\lambda$','interpreter','latex','FontSize',14)
% -------------------------------------------------------------------------
% Ende Programm
% -------------------------------------------------------------------------

