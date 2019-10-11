
function [info] = fast_read_CMN3(year,DOY,Time,d_Time, Elevation,latstart,latend,lonstart,lonend,numbstatn,satheight,min_dist_stn,Data_location)

% read the CMN files and extract the reciever position 
% the reciever coordinates are located on line three 
% output in  altitude,lat,lon 
% Also extract all the needed parameters like STEC VTEC Elmat and azimuth
% ... 
% By Nicholas Ssessanga while at Chungnam National University
% 2019-july-02

% this function reads data from a whitespace delimited ascii file
% with colum headings
% it tries to work out whether the columns contain numbers or text   
% check which system we are running on 
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
    LghtCMN = length(dir([Data_location,slashx,'*',DOY_folder,'-',Year_fld,'*']));
    % get only data corresposing to the  day and year we are processing 
    CMNfile = dir([Data_location,slashx,'*',DOY_folder,'-',Year_fld,'*']);
	
    
    % this is the range of the needed time in decimal hours 
    ndtime =Time;
    ndtime1 = ndtime -(60*d_Time)/(3600);
    ndtime2 = ndtime +(60*d_Time)/(3600);
 
    I = 0;
    inuse = 1; 
    K = 1;
    Lnext=0;
    p=1;
    
     while 1
      
            I = I+1;
            % check and see if maxmum number of stations is reached 
            % if so, break and exit 
            if (I > LghtCMN)||( K > numbstatn)
                break
            end 
            file_name = [Data_location,slashx,CMNfile(I).name];
           % fprintf('Reading "%s"\n', file_name)
            
            % open a data file and look for the required data 
            [fid, message] = fopen(file_name);
            
            if (fid == -1)
               
                error(message);
            end
             
            % get the Rx position altitude above sea level in meters,
            % Lat, Log  in degerees
            % located on the 3rd line 
            m = 1;
            while 1 % while loop to get receiver position 
                line = fgetl(fid);
            
                if ~ischar(line)
                    Rx = [ ];
                    break
                end
                if m>3, break, end
                if (m==3)
                    [t] = sscanf(line, '%f %f %f');
                    Rx = [t(3) t(1) t(2)];

                end
                m=m+1;
           end % end while loop to get receiver position 
           
           % reciever position is not available do nothing  
           if  ~ isempty (Rx)
                rlat = Rx(2);
                rlon = Rx(3);
             
              % check if this station is within the region of interest
              if  ( rlat>= latstart) && (rlat<= latend)&& (rlon >= lonstart) && (rlon <= lonend)

                   % check if this station is outside the minimum distance 
                    if K > 1
                        [inuse]=checkmindistance ( Rx, Rxo,min_dist_stn);
                    end 



                    % if this station is outside the minimum distance continue 
                    % processing
                    % Also if K = 1, then this is the intial station so we do not
                    % need to satisfy the minimum distance requirement. 

                    if inuse==1 || K ==1 
                        header_line = fgetl(fid); % reads line without end of line character
                        % get the Position in open file (in bites)
                        data_start_position = ftell(fid);
                        % get the labels for each column
                        % names = regexp(header_line, '\S*', 'match');
                        % now try parsing the first data line
                        % this will be used as dummy to set up the data
                        % format 
                        test_line = fgetl(fid); 
                        % split up the data line into strings 
                        test_data = regexp(test_line, '\S*', 'match');
                        % This should give us the number of data columns 
                        num_cols = length(test_data);
                        format_str = '';
                        all_numbers = 1;
                        
                        % now create the format of the data using regexp
                        % We use the first line as an example to formulate the data format 
                        % here every thing is a float %f 
                        for i = 1: num_cols
                            match = regexp(test_data{i}, '([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?', 'match');
                            if isempty(match)
                                format_str = [format_str '%s'];
                                all_numbers = 0;
                            else
                                format_str = [format_str '%f']; % here %f %f %f %f %f %f %f %f %f %f
                            end
                        end 
                        % after creating the data format, now lets read the
                        % whole data in the file. 
                        % Move to specified position we were before in file using the fseek

                        fseek(fid, data_start_position, 'bof');
                        % choose whether to use textscan or fread 
                        % if the data contains strings then use text scan
                        % otherwise use fread for floating data  
                        if all_numbers == 0
                            data = textscan(fid, format_str);
              
                        else
                
                          %  raw_data = fread(fid,'uint8=>char');
                          % data_array = reshape(sscanf(raw_data, format_str), num_cols, []);
                            data = textscan (fid, '%*f %f %f %f %f %*f %*f %f %f %*f', 'HeaderLines',5);
                            fclose(fid);
%                             data = {};
%                             for i = 1: num_cols
%                                 data{i} = data_array(i, :)';
%                              
%                             end
                            
                            %JD = out{1}; 
                            Time1 =data{1};%data_array(2, :);% data{2};%2; 
                            PRN   =data{2};%data_array(3, :); %data{3};%3; 
                            Azmat  =data{3};%data_array(4, :);% data{4};%4; 
                            Elmat  =data{4};%data_array(5, :);% data{5};%5;  
                            %Lat= out{6}; 
                            %Long = out{7}; 

                            STEC =data{5};%data_array(8, :);% data{8}; %8; 
                            VTEC=data{6}; %data_array(9, :); %data{9}; %9;
                            
                            % free up memory
                            clear data

                            % Now extract the required Time with Elmatation angle greater than minmum to
                            % avoid any multipath, also avoide any bad STEC that
                            % might be negative
                            Time1index = ((roundn(ndtime1,-5) <= Time1) & (Time1 <= roundn(ndtime2,-5))...
                                                    &(abs(Elmat) ~= 90)&(abs(Elmat) > Elevation) & (STEC > 0)& (PRN==15 | PRN==20));

                            if (~isempty(PRN(Time1index))&& ~isempty(Time1(Time1index))&&  ~isempty(STEC(Time1index))...
                                                            && ~isempty(Elmat(Time1index))&&  ~isempty(Azmat(Time1index)))



                                    PRN = PRN(Time1index);    % get corresponding PRN
                                    Time1 = Time1(Time1index); % get corresponding Time
                                    STEC = STEC(Time1index);  % get corresponding stec 
                                    Elmat = Elmat(Time1index);     % get corresponding Elmatation matrix
                                    Azmat =Azmat(Time1index);     % get corresponding azimuth matrix
                                    VTEC=VTEC(Time1index);  
                                    stationLat =  Rx(2);          % get corresponding receiver coordinates(Lat) 
                                    stationLong = Rx(3);          % get corresponding receiver coordinates (Long)
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
                                             info.RxLat(p:Lnext)  = B(:,2);
                                             info.RxLon(p:Lnext)  = B(:,3);
                                             info.RxAlt(p:Lnext)  = B(:,1);
                                             info.prn(p:Lnext)    = PRN(ind);
                                             info.time(p:Lnext)   = Time1(ind);
                                             info.vtec(p:Lnext)   = VTEC(ind);  

                                             p =  Lnext+1; 
                                             
                                             Rxo =Rx; % set the previous station to the new station
                                             K=K+1;   % increase the station counter 
                                             inuse = -1; % the the inuse flag
                                     end % end for isempty ind

                                    

                            end
                            
                        end

                    end 


              else 
                  fclose (fid);
                  
                
                  % put here code for not in the region of interest  
                  
              end  % end of if in selected region of interest



           else % this file has no receiver position  
                % close the fid and move to the next file  
             fclose (fid);  
          end 
          
              
     end  % end of while loop for all files

end 

     
function [inuse]=checkmindistance ( Rx, Rxo,min_dist_stn)
    % This function checkes if the next station to be analysed is to close to the previously 
    % analysed station. To avoide redundancy, we neglect any station that i s located within
    % the minimum distance area, currently set to 100 km ~ 1 deg 
    dist_deg  = ((Rx(2)-Rxo(2))^2 + (Rx(3)-Rxo(3))^2)^.5;

    if dist_deg < min_dist_stn 
        inuse= 0 ;
    else 
         inuse= 1; 
    end 

end 
     
     
