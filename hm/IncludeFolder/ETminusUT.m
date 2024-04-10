% -------------------------------------------------------------------------
% ETminusUT.m
% -------------------------------------------------------------------------
% MATLAB-Programm zum Kapitel Himmelsmechanik aus
% "Physikalische Fingerübungen" von Michael Kaschke und Holger Cartarius
% unter Mitwirkung von Ulrich Potthoff
% Alle Rechte bei den Autoren
% Freier Gebrauch mit Buch und/oder Angabe der Quelle erlaubt.
% -------------------------------------------------------------------------
% Berechnet Delta T = ET - UT als Funktion des Julianischen Datums 
% über Polynome (anwendbar nur im Zeitraum 1960 -2150)
% Berechnung von DeltaT über Polynome nach Espenak
% https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
% -------------------------------------------------------------------------
% Input:  T Julianisches Datum
% Output: DeltaT in sec
% -------------------------------------------------------------------------

function DelTSec = ETminusUT(T)
  dt = datetime(T,'ConvertFrom','juliandate');
  y  = dt.Year;
  if y < 1961 % vor 1960
        DelTSec = 30;
  else 
    if y < 1985 % Between years 1961 and 1986 calculate where: t = y - 1975
        t = y - 1975;
        DelTSec = 45.45 + 1.067*t - t^2/260 - t^3 / 718;
    else
      if y < 2005 % Between 1986 and 2005 calculate where: t = y - 2000
        y = y -2000;
        DelTSec = 63.86 + 0.3345*y - 0.060374 *y^2 + 0.0017275*y^3 + ...
                      0.000651814*y^4 + 0.00002373599*y^5;
      else
        if y < 2051 % Between years 2005 and 2050 calculate 
                    % where: t = y - 2000
            y = y - 2000;
            DelTSec = 62.92 + 0.32217*y + 0.005589*y^2;
        else
          if y < 2150 % Between years 2050 and 2150, calculate:
            DelTSec = -20 + 32 * ((y-1820)/100)^2 - 0.5628 * (2150 - y);
          else 
            DelTSec = 220;
          end
        end
      end
    end 
  end
end

% -------------------------------------------------------------------------
% Ende Funktion
% -------------------------------------------------------------------------