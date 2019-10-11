function  [hcorrPtr]= ComputeHorizontalCorr (All_lat,All_HGT,All_lon)
  
%*******************************************************
% This function computes the correlations between the points
% in the horizontal grid plane.  Pay attention to the units
% in self.latCorrUnits and self.lonCorrUnits which determine
% how we handle distances.
%********************************************************


useLatCorr = 1;
useLonCorr = 1;
%useAltCorr = 1;

% Some tolerances
ExpMax = 20.0;
small = 0.0000001;


[magcord] = magnetvec(All_lat,All_HGT,All_lon);
numPlane = length(magcord(:,1));
gammamin = 10000000.00; 
% solar terminator to be revised later 
% to include whether a point is in night or sun side 
% for the start all points lie in the same zone

dailyIndexVecPtr = ones(1,numPlane);
hcorrPtr = zeros(numPlane,numPlane);

    for i=1 :numPlane
       magLat1 =  magcord (i,1);
       magLon1 = magcord (i,2);
      % CALL LatLonToRangeAzi(magCoordlat,magCoordLon,magLat1,magLon1,gammavec,azvec)
       [gammavec,~] = LatLonToRangeAzi_vector(magcord(:,1), magcord(:,2), magLat1,magLon1);
       gamma = min(gammavec(gammavec > 0.0001));
       if (gamma < gammamin) 
           gammamin = gamma;
       end 
    end
    
    
    if (gammamin < 2.0*pi/180)
        gammamin =  2.0*pi/180;
    end

% Minimum correlation value allowed - used to default to in certain cases
    L2_small = (gammamin)^2.0;
   for i=1:numPlane 
       magLat1 =  magcord (i,1);
       magLon1 = magcord (i,2);
      l_scale_1 = 5*pi/180;%GetLatScale(this, magLat1, i)
      lon_scale_1 = 8*pi/180;%GetLonScale(this, magLat1, i)   
      
      
      for j=1:numPlane
       xx = 0.D0;
       magLat2 =  magcord (j,1);
       magLon2 = magcord (j,2);
         l_scale_2 = 5*pi/180;%(this, magLat2, j)
         lon_scale_2 =8*pi/180;% GetLonScale(this, magLat2, j)
         % get the azimuths from 1-2 and from 2-1. great circle distance is the same but
         % azimuths are not 180 degrees apart on a sphere% so for symmetry need both.
         [~,az1] = LatLonToRangeAzi_vector( magLat1, magLon1 , magLat2,  magLon2);
         [gamma,az2] = LatLonToRangeAzi_vector(magLat2,  magLon2, magLat1, magLon1);
        
         
         % for the case of both lat and lon correlations, we have to 
         % adjust each correlation length by the direction.
         oneByL2 = (((cos(az1))^2)/l_scale_1^2) + ((sin(az1))^2)/(lon_scale_1^2);
         Ll1 = sqrt(1.0/oneByL2);
         oneByL2 = (((cos(az2))^2)/l_scale_2^2) + ((sin(az2))^2)/(lon_scale_2^2);
         Ll2 = sqrt(1.0/oneByL2);

         % for the case where we only have lats or lons
         Llat2 = l_scale_1*l_scale_2;
         Llon2 = lon_scale_1*lon_scale_2;
        % crossing terminator - small
         if (dailyIndexVecPtr(i) ~= dailyIndexVecPtr(j)) 
            Llat2 = L2_small;
            Llon2 = L2_small;
            Ll1 = sqrt(L2_small);
            Ll2 = sqrt(L2_small);
         end
  
 
% slect which correlations to compute (lon, lat or alt)
         if ((useLatCorr==0)&&(useLonCorr==1)) 
            if(abs(magLat1 - magLat2)< small)
               xx = ((gamma)^2)/Llon2;
            end 
         end 
         if ((useLatCorr==1)&&(useLonCorr==0)) 
            if (abs(magLon1 - magLon2)<small) 
               xx = ((gamma)^2)/Llat2;
            end
         end
         if((useLatCorr==1)&&(useLonCorr==1)) 

            xx = (gamma)^2/(Ll1*Ll2);
         end 
         if (xx>ExpMax)
             xx = ExpMax;
         end 
        hcorrPtr(i,j) = exp(-1.*xx);

      end 
% Check the diagonal terms and reset just in case
      hcorrPtr(i,i) = 1.0  ;    
   end 
