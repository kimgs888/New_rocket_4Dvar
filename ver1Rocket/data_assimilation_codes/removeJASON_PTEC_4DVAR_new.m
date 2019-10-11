function [real_STEC] =  removeJASON_PTEC_4DVAR_new(Elv,STEC,VTEC,sun_up_down)
% Our reconstruction altitude is capped 100- 1336km; This is the
   % the JASON-1 Altitude. However GPS orbit altitude is ~ 20000km, 
   % therefore STEC includes plasmasphere(plm). 
   % to remove the contribution of plm to STEC before assimilation
   % we use studies in a paper by  E Yizengaw - ?2008  using data from
   % JASON. at mid latitude 30% and  15 % of VTEC during
   % night and day  respectively
  % The time from file is in Julian so covert it to normal time 
  % for easy reading and management 
 %[ day month year hour minu sec ] = Julian2Greg(time1) ;

   % hourX = (hour + minu/60+ sec/3600);
  % assume night time starts at 10 UT in the northern Hemisphere (Korea are)
 
    nighttimeidx = (sun_up_down<=0);
    daytimeidx = (sun_up_down>=1);

    thetanight = sin(Elv(nighttimeidx)*pi/180);

    
    real_STEC = zeros(size(STEC));
    real_STEC(nighttimeidx) = STEC(nighttimeidx)-(0.3*VTEC(nighttimeidx)./thetanight);
    % instead of using 15% for day I  have used 10% as an average from
    % values given in the paper, I think we are over estimating this values
    % if we use 15 % according to the results 
    thetaday= sin(Elv(daytimeidx)*pi/180); 
   
    real_STEC(daytimeidx) = STEC(daytimeidx)-(0.15*VTEC(daytimeidx)./thetaday);