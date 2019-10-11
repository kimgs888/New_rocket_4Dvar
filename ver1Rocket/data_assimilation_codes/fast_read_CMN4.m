
function [info] = fast_read_CMN4(year,DOY,Time,d_Time, Elevation,latstart,latend,lonstart,lonend,numbstatn,satheight,min_dist_stn,Data_location)

% read the CMN files and extract the reciever position 
% the reciever coordinates are located on line three 
% output in  altitude,lat,lon 
% Also extract all the needed parameters like STEC VTEC Elmat and azimuth
% ... 
% By Nicholas Ssessanga while at Chungnam National University
% 2019-july-02


    if ispc 
             slashx = '\';
        elseif isunix 
            slashx = '/';
    end 
 
  
    % fomulate the folder  and path to the data   
     Year_fld = num2str(year); 

    if DOY < 10 

        DOY_folder = ['00',num2str(DOY)];

    elseif   DOY >=10 && DOY <100

        DOY_folder = ['0',num2str(DOY)];

    elseif DOY >= 100

         DOY_folder = num2str(DOY);

    end 
    
    % Get the number of files to be processed 
    % remember each file is one station
    
    
    % Fomulate the commom file string in all files to be processes
  
    Common_str = ['*',Year_fld,DOY_folder,'.bin'];

    % Get all files with that common string 
    datafiles = dir([Data_location,slashx,Common_str]);

    
    % Total number of files selected  
    Totalfiles = length(datafiles);
    

    
    % Center the required time around the length of the collection window 
    % within this period all the ionosphere is assumed to be stationary 
    % \...ndtime.../
    ndtime =Time;
    ndtime1 = ndtime -(60*d_Time)/(3600);
    ndtime2 = ndtime +(60*d_Time)/(3600);
 
    % intialize all parameters 
    openedfiles = 1; % opened files counter 
    
    inuse = 1; 
    numbstatni= 1;
    Lnext=0;
    p=1;

   % check and see if maxmum number of stations is reached 
   % if so, break and exit 
while ((openedfiles <= Totalfiles)&&( numbstatni <= numbstatn))
      
        

        file_namei = [Data_location,slashx, datafiles(openedfiles).name];
        
        % check and see if this file does exist in the data folder  
        % if it does, "exist" outs a int 2 
        
        if exist(file_namei,'file')== 2
            
                % open that data file and look for the required data 
                % it a read only, so include 'r'
                [fid, ~] = fopen(file_namei,'r');

                if (fid == -1) % if the file identifier (fid) is -ve that file cant be opened 

                         disp([file_namei, ' is corrupted or does not exist'])
                
                
                
                
                else % otherwise all is fine, so lets proceed
                    
                       % The first 3 bytes in this binary file 
                       % describe the receiver location  Lat Lon and
                       % altitude 
                        Rx    = fread(fid,[1,3],'*double');
                     

                       % if reciever position is not available, we can not
                       % use this file, so close everthing and move on the
                       % next file 
                     
                       if  ~ isempty (Rx)
                            rlat = Rx(1);
                            rlon = Rx(2);

                          % if Rx is available check if this station is within the region of interest
                          if  ( rlat>= latstart) && (rlat<= latend)&& (rlon >= lonstart) && (rlon <= lonend)

                               % if Rx  is in region of interest check if its not to close 
                               % to the previously analysed station. 
                               % However only run this check if we have
                               % processed the first station
                               
                                if  numbstatni > 1
                                    [inuse]=checkmindistance ( Rx, Rxo,min_dist_stn);
                                end 



                                % if this station is outside the minimum distance continue 
                                % processing
                                
                                % Also if numbstatni = 1, then this is the intial station so we do not
                                % need to satisfy the minimum distance requirement. 

                                if inuse==1 || numbstatni ==1 
                                
                                        
                                   % Read the binary data  
                                   data_array = fread(fid,'*double');
                                   fclose(fid);
                                      
                                   % reshape the data into the number of
                                   % data columns in cmn files these are
                                   % 10 column
                                   % [], n divide the data into n equal
                                   % columns 
                                   
                                   data_array = reshape(data_array,[],6);

                                   Time1 = data_array(:,1);%data_array(:,2);%data{1};%data_array(2, :);% data{2};%2; 
                                   PRN   = data_array(:,2);%data_array(:,3);%data{2};%data_array(3, :); %data{3};%3; 
                                   Azmat = data_array(:,3);%data_array(:,4); %data{3};%data_array(4, :);% data{4};%4; 
                                   Elmat =data_array(:,4);% data_array(:,5);%data{4};%data_array(5, :);% data{5};%5;  
                                   %Lat= out{6}; 
                                   %Long = out{7}; 
                                   STEC = data_array(:,5);%data_array(:,8);%data{5};%data_array(8, :);% data{8}; %8; 
                                   VTEC = data_array(:,6);%data_array(:,9);%data{6}; %data_array(9, :); %data{9}; %9;
                                   % free up memory
                                   clear data data_array

                                        % Now extract the required Time with Elmatation angle greater than minmum to
                                        % avoid any multipath, also avoide any bad STEC that
                                        % might be negative
                                        Time1index = ((roundn(ndtime1,-5) <= Time1) & (Time1 <= roundn(ndtime2,-5))&(abs(Elmat) ~= 90)&(abs(Elmat) > Elevation) & (STEC > 0));

                                        if (~isempty(PRN(Time1index))&& ~isempty(Time1(Time1index))&&  ~isempty(STEC(Time1index))...
                                                                        && ~isempty(Elmat(Time1index))&&  ~isempty(Azmat(Time1index)))



                                                PRN = PRN(Time1index);    % get corresponding PRN
                                                Time1 = Time1(Time1index); % get corresponding Time
                                                STEC = STEC(Time1index);  % get corresponding stec 
                                                Elmat = Elmat(Time1index);     % get corresponding Elmatation matrix
                                                Azmat =Azmat(Time1index);     % get corresponding azimuth matrix
                                                VTEC=VTEC(Time1index);  
                                                stationLat =  Rx(1);          % get corresponding receiver coordinates(Lat) 
                                                stationLong = Rx(2);          % get corresponding receiver coordinates (Long)
                                                satalitude =satheight;        %  assume a transmitter (satellite) at this height


                                                % call calcsatposnk  to calculate the satellite position (Lat, Long) using Elmatation and azimuth 


                                                [sat]= calcsatposnk (Elmat, Azmat, stationLong,stationLat,satalitude);

                                                % Now limit the outcome to rays within the grid 

                                                ind = ((sat.lat>=latstart) & (sat.lat<=latend)& (sat.long >=lonstart) & (sat.long <=lonend));

                                                if ~isempty(sat.lat(ind)) %check to see if index is empty if so no rays in that file 

                                                         B = repmat(Rx,size(sat.lat(ind)));
                                                         C = repmat(satalitude,size(sat.lat(ind))); 

                                                         Lnext = Lnext + length(sat.lat(ind)); 

                                                         info.satlat(p:Lnext) = sat.lat(ind);
                                                         info.satlon(p:Lnext) =  sat.long(ind);
                                                         info.Satalt(p:Lnext) = C(:);
                                                         info.elv(p:Lnext)    = Elmat(ind);
                                                         info.Azm(p:Lnext)    = Azmat(ind);
                                                         info.stec(p:Lnext)   =STEC(ind);
                                                         info.RxLat(p:Lnext)  = B(:,1);
                                                         info.RxLon(p:Lnext)  = B(:,2);
                                                         info.RxAlt(p:Lnext)  = B(:,3);
                                                         info.prn(p:Lnext)    = PRN(ind);
                                                         info.time(p:Lnext)   = Time1(ind);
                                                         info.vtec(p:Lnext)   = VTEC(ind);  

                                                         p =  Lnext+1; 

                                                         Rxo =Rx; % set the previous station to the new station
                                                         numbstatni=numbstatni+1;   % increase the station counter 
                                                         inuse = -1; % the the inuse flag
                                                 end % end for isempty ind



                                        end

                                    

                                end 


                          else  % this station is not in the region of interest 
                                % close file and move on to the next
                                
                               fclose (fid); 


                              

                          end  % end of if in selected region of interest



                       else % this file has no receiver position  
                            % close the fid and move to the next file  
                         fclose (fid);  
                      end 

               
                end % end if data file can be opened 
               
        else  % this binary file does not exist 
            disp([file_namei, ' does not exist!'] )

        end % check if .bin file exist in data folder  
         
         
      % increment number of opened files 
      openedfiles = openedfiles+1;
end  % end of while loop for all files

end 

     
function [inuse]=checkmindistance ( Rx, Rxo,min_dist_stn)
    % This function checkes if the next station to be analysed is to close to the previously 
    % analysed station. To avoide redundancy, we neglect any station that i s located within
    % the minimum distance area, currently set to 100 km ~ 1 deg 
    dist_deg  = ((Rx(1)-Rxo(1))^2 + (Rx(2)-Rxo(2))^2)^.5;

    if dist_deg < min_dist_stn 
        inuse= 0 ;
    else 
         inuse= 1; 
    end 

end 
     
     
