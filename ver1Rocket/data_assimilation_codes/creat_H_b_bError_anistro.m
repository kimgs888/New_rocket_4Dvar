function [dH,db,bError, Rxinfo]=creat_H_b_bError_anistro (year,DOY,Time,d_Time,Elevation,numbstatn,satheight,...
                                                                                    All_HGT,All_lat,All_lon,min_dist_stn,Data_location)
                                               

% This functions extracts delta geometry(dH) and delta stec(db)
% in an observation arc, we subtract the minmum stec from all other rays 
% In this process we dont need to calculate the reciever and satellite
% biases
% However, since we are using absolute STEC, the subtruction process is neglected

% By Nicholas Ssessanga while at Chungnam National University
% 2019-july-02

    % Get all the data in the selected range 
    RE = 6371e3; 
    latstart = min(All_lat);
    latend = max(All_lat);
    lonstart= min(All_lon);
    lonend = max(All_lon);
                            % ~~~~~~  Extract the data structure from the Data files  ~~~~~~ 
                            % ~~~~~~  the all info needed from the data files is recored and stored in data structure infm 
    [infm]= fast_read_CMN5(year,DOY,Time,d_Time, Elevation,latstart,latend,lonstart,lonend,numbstatn,satheight,min_dist_stn,Data_location);
                                   
  
 
    stlat = infm.satlat;   % get the satellite latitudes 
    stlon = infm.satlon;   % get the satellite longitudes 
    stalt = infm.Satalt;   % get the satellite altitudes
    Rxlat = infm.RxLat;    % get the Receiver  latitude
    Rxlon = infm.RxLon;    % get the Receiver  Longitude
    Rxalt = infm.RxAlt;    % get the Receiver  altitudes
    % remove these fields since we are done with them 
    infm = rmfield(infm, 'satlat');
    infm = rmfield(infm, 'satlon');
    infm = rmfield(infm, 'Satalt');
    infm = rmfield(infm, 'RxLat');
    infm = rmfield(infm, 'RxLon');
    infm = rmfield(infm, 'RxAlt');

    %     ~~~~~~~ perform the ray tracing ~~~~~~~ 

    Rx = sphcart([RE+Rxalt',[Rxlat',Rxlon']*pi/180]);
    Tx = sphcart([(RE+stalt)',[stlat',stlon']*pi/180]);

    GHT = (All_HGT.*1E3)+RE;
    Grid  = setgrid('Rad',GHT,'Lat',All_lat,'Lon',All_lon);
    Ray   = gpsray(Rx,Tx,Grid);
        
    % free up memory space and  and do ray tracing 
    clear Rx Tx stlat stlon stalt 
    [H] = raytracer(Grid, Ray);
    % clear variables used in ray tracing 
    clear  Grid Ray 
        
    %   ~~~~~  end of ray tracing  ~~~~~~~ 
    
    %    ~~~~~  get corresponding variables ~~~~~~~ 
        prn = infm.prn';
        TIME = infm.time';
        vtec = infm.vtec';
        elev = infm.elv';
        STEC = infm.stec';
    
    %    ~~~~~~~ Remove PTEC ~~~~~~~ 
    
    % Use this sun_up_down to subtract the PTEC 
    % sun_up_down determines weather the station is in day or night region 
    % then deduct the plasmasphere accordingly 
     [sun_up_down] = Day_or_night_fxn (Time,DOY,Rxlat',Rxlon', Rxalt'); 
%  
%     [STEC] =  removeJASON_PTEC_4DVAR_new(elev,STEC,vtec,sun_up_down);
    b = STEC*1E16; % data without PTEC

    clear infm real_STEC sun_up_down
    %    ~~~~~~~ Done Removing  PTEC ~~~~~~~ 
    
    %    ~~~~~~~ remove all the bad Rays ~~~~~~~ 
    
    % if the total length of the ray is less than 1e-17
    % or STEC is less than .5 TEC 
    % thats a bad ray.
    
%     % H bad 
    F = sum(H,2)
%     indF = (F<=1e-17);
%     clear F 
%     if any(indF)
%         
%         H(indF,:)=[];
%         b(indF)  =[];
%         Rxlat(indF)=[];   
%         Rxlon(indF)=[];  
%         prn(indF)=[];
%         TIME (indF)=[];
%         clear inF
%     end 
%     % b bad 
    Fb =b/1e16;
    indFb = (Fb<=0.5 | Fb>=500);
    clear Fb 
    if any(indFb) 
  
        H(indFb,:)=[];
        b(indFb)=[];
        Rxlat(indFb)=[];   
        Rxlon(indFb)=[];   
        prn(indFb)=[];
        TIME (indFb)=[];

        clear inFb
    end
    %    ~~~~~~~ Done removing all the bad Rays ~~~~~~~ 

    
%% we have all the info we need, now start the dH db process 



    srtrng = 1; 
    endrng =0;
    number_base_rays = 0;
    %rcvno =1 ;
    while 1 
    if isempty(TIME)

        dH= dH'; 
        db= db'; 
        bError = bError';
        break

     else

        Rxloni = Rxlon(1);
        Rxlati =  Rxlat(1);
        prni = prn(1);
        %flagi=flag(1);

        drxlon = abs(Rxlon'-Rxloni);
        drxlat = abs(Rxlat'-Rxlati);
        dprn   = abs(prn-prni);
        %dflag  = abs(flag-flagi);
        indx = ((drxlon<0.00001)& (drxlat<0.00001)&(dprn<0.00001));%&(dflag<0.00001));

        %  '############################################'   
        %   [Rxlon(indx) Rxlat(indx) prn(indx) flag(indx) TIME(indx) b(indx)/1e16]
        %'******************************************' 

        b_slect = b(indx);


        if size(b_slect,1)>0
            number_base_rays =number_base_rays +1;
            H_slect = H(indx,:)';
            b_error = b_slect-mean(b_slect);

            endrng = endrng + size(H_slect,2);
            dH(:,srtrng:endrng)= H_slect;
            db(:,srtrng:endrng)= b_slect'; 
            bError(:,srtrng:endrng)= b_error'; 
            Rxinfo.RxLat(srtrng:endrng) =Rxlat(indx);
            Rxinfo.RxLon(srtrng:endrng)= Rxlon(indx) ;
            Rxinfo.prn(srtrng:endrng)=prn(indx);
            Rxinfo.time(srtrng:endrng)=TIME(indx);

             srtrng = endrng+1; 

            Rxlon(indx) = [];
            Rxlat(indx) = [];
            prn(indx) = [];
            %flag(indx) = [];
            TIME(indx) = [];
            b(indx)=[];
            H(indx,:)=[];
            %db = [db;dbi];
            %dH = [dH;dHi];
          else
            Rxlon(indx) = [];
            Rxlat(indx) = [];
            prn(indx) = [];
            %flag(indx) = [];
            TIME(indx) = [];
            b(indx)=[];
            H(indx,:)=[];
        end 
    end 
    end 
end 



