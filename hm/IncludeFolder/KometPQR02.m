% -------------------------------------------------------------------------
% KometPQR02.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Optimierung Kometenbahn und PQR
% -------------------------------------------------------------------------

function KometData=KometPQR02(KometDataIn, T, BaPaK,HK)
  KometData = KometDataIn;
  for n=1:length(KometData.M)
     if (abs(KometData.M(n)) > 0.2) 
           KometData.xyz(:,n)=KometData.xout_ell(:,n);
           KometData.v(:,n)=KometData.vout_ell(:,n);
      else
           KometData.xyz(:,n)=KometData.xout_par(:,n);
           KometData.v(:,n)=KometData.vout_par(:,n);
     end
  end
  OmegaPn=deg2rad(BaPaK.OmegaP(HK));
  omegaPn=deg2rad(BaPaK.omegaP(HK));
  iPn=deg2rad(BaPaK.iP(HK));
  % PQR Multiplikation
  yin=[KometData.xyz(1,:);KometData.xyz(2,:);KometData.xyz(3,:)];
  yout=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
  KometData.xyz=yout;
  % Geschwindigkeitsvektoren heliozentrisch des Kometen
  yin=[KometData.v(1);KometData.v(2);KometData.v(3)];
  yout=mtimes(PQR(-OmegaPn,-iPn,-omegaPn),yin);
  KometData.v=yout;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------