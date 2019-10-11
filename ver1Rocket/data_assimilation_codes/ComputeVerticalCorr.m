function [vCorrPtr ]= ComputeVerticalCorr(All_lat,All_HGT,All_lon)
  
%***********************************************************
% This function computes the vertical correlations between altitudes
% for the standard grid vector on which everything is built.  It is assumed
% that this vector is the same at all sites.  All values have previously been interpolated
% to this grid.  The input htvec is the vector of altitude grid values.
%
% Pay attention to the length scale units in self.lengthScaleUnits - these determine if we
% are using absolute indices or distances.
% *********************************************************
% correlation length    


% Set up some initial tolerance
 small = 0.01;
 ExpMax = 20.0;

% How many loops in heights?
%[magcord] = magnetvec(All_lat,All_HGT,All_lon);

htVec = All_HGT;
htVec = htVec(:);
numHeights = length(htVec);
% If the units are indices, build an array of indices to represent the heights,
% otherwise get the actual heights.
correl = [60:100,405:1000,1050-2000,2200:3000,4000:20000];
scaleLegth=[ones(size(60:100))*100, ones(size(405:1000))*200,...
           ones(size(1050-2000))*400,ones(size(2200:3000))*750,...
           ones(size(4000:20000))*1500];
vCorrPtr = [];%zeros(numHeights,numHeights);
hmax = 450;


          for i =1:numHeights 
             i_hght = htVec(i);
             [~,indexi]  = min(abs(correl(:)- i_hght) );
            for j=1:numHeights
              j_hght = htVec(j);
              [~,indexj]  = min(abs(correl(:)- j_hght) );
              lht2 = scaleLegth(indexi)*scaleLegth(indexj);
              if(lht2<small) 
                  lht2 = small;
              end
              dx2 = (htVec(j) - htVec(i))^2;
              xx = dx2/lht2;
              if (xx>ExpMax) 
                  xx = ExpMax;
              end
              if ((htVec(j) > hmax && htVec(i) < hmax) || ...
                  (htVec(i) > hmax && htVec(j) < hmax))
               vCorrPtr(j,i) = -1.*exp(-1.0*xx);  
              else
                vCorrPtr(j,i) = exp(-1.*xx);
              end
            end
    
         end

