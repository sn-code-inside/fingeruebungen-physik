% -------------------------------------------------------------------------
% SAchse.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Fingerübungen der Physik" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnung der Schattenachse der SoFi
% -------------------------------------------------------------------------

function SA = SAchse(MoonKoord, SunKoord)
   fac = 0.996633; % Abplattung Erde
   DM  = 0.5450;   % Monddurchmesser in Erdradien
   DS  = 218.25;   % Sonnendurchmesser in Erdradien
   X_M = MoonKoord.equ1; % kartes. äquatoriale Koord. v. Mond und Sonne
   X_S = SunKoord.equ1;
   X_M(3) = X_M(3)/fac; % Streckung der z-Koordinate wg Abplattung
   X_S(3) = X_S(3)/fac; 
   X_MS   = X_M-X_S;
   R_MS   = norm(X_MS);
   E_XS   = X_MS/R_MS; % Richtungsvektor Sonne-Mond
   p0     = -sum(X_M .* E_XS); % Entfernung Mond-Hauptebene
   Dqua   = p0*p0 + 1.0-norm(X_M)^2;
   r0     = sqrt(1-Dqua); % Entfernung Erdmitte-Schattenachse
   DKS    = (DS-DM)*p0/R_MS-DM; % Durchmesser Kernschatten auf Hauptebene
   DHS    = (DS+DM)*p0/R_MS+DM; % Durchmesser Halbschatten auf Hauptebene
   XB = [NaN;NaN;NaN];
   if r0 < 1.0 % Schattenachse trifft auf Erde
      p   = p0 -sqrt(Dqua);
      DKS = (DS-DM)*(p/R_MS)-DM; % Durchmesser Kernschatten auf Erdoberfläche
      XB  = X_M + p*E_XS;        % kartesisches Koordinaten des Schattenkegels
      XB(3) = XB(3)*fac;         % Rückrechnung der Abplattung
      if DKS > 0 
          Phase = "ringförmig";
      else
          Phase = "total";
      end
   else
      if r0 < (1.0+0.5*abs(DKS))
         if DKS > 0 
             Phase = "ringförmig, nicht zentral";
         else 
             Phase = "total, nicht zentral";
         end
      else
          if r0 < (1.0+0.5*abs(DHS))
             Phase = "partiell";
          else
             Phase = "keine Finsternis";
          end
      end 
   end
   SA.Phase = Phase;
   SA.xyz   = XB;
   SA.equ   = CalcAnglesfromXYZ(SA.xyz);
   SA.E_XS  = E_XS;
   SA.DKS   = DKS;
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------
