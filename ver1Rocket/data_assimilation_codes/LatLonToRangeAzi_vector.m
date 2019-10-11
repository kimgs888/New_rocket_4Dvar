function [range2, az] = LatLonToRangeAzi_vector(lati, loni, lato, lono)
      
    



      range2 = zeros(length(lati));
      az =zeros(length(lati));

      npts = length(lati);

      ok =  (length(lati) == npts )& (length(loni) == npts) & (length(range2) == npts) & (length(az) == npts);

      if (ok ~= 1)
      disp('Ray:LatLonToRangeAzi_vector(): The input lat,lon, range, and az vectors are not the same length')
      else

       pib2 = pi/2;
       npts =length(lati);

        if (npts > 0) 

          for i = 1:npts

         % n_pole test case

            if (abs(lati(i) - pib2) < 1e-6)

              range2(i) = pib2 - lato;
              az(i) = pi;

            % s_pole test case

             elseif (abs(lati(i) + pib2) < 1.0D-6) 

              range2(i) = pib2 + lato;
              az(i) = 0.0D0;

            % the same

             elseif ((abs(lati(i) - lato) < 1.0D-6) &&...
                    ( abs(loni(i) - lono) < 1.0D-6)) 

              range2(i) = 0.0D0;
              az(i) = 0.0D0;

            % general

            elseif (abs(lati(i) - pib2) > 1.0D-6  &&...
                     abs(lati(i) + pib2) > 1.0D-6  &&...
                     (abs(lati(i) - lato) > 1.0D-6||...
                      abs(loni(i) - lono) > 1.0D-6))
              sin1 = sin(lati(i));
              cos1 = cos(lati(i));
              s2 = sin(lato);
              cr = sin1 * s2 + cos1 * cos(lato) * cos(lono - loni(i));

              range2(i) = acos(cr);

              ca = (s2 - sin1 * cr) / (cos1 * sin(range2(i)));
              if (ca < -1.0D0) 
                  ca = -1.0D0;
              end
              if (ca >  1.0D0) 
                  ca =  1.0D0;
              end

              az(i) = acos(ca);

              if (sin(lono - loni(i)) < 0.0D0) 
                az(i) = az(i) - 2.0D0 * pi;
              end
             end
         end

        end
     end


