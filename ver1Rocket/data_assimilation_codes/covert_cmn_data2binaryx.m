function [] = covert_cmn_data2binaryx ( Data_location,year,doy,All_lat,All_lon,numbstatn,min_dist_stn)    
% This code coverts the *.CMN files into Binary format 
% Its way quicker to read binary data than Ascii data format 
% The computation time is highly reduced when we deal with binary files 
% By Nicholas Ssessanga while at CNU 

latstart = min (All_lat);
latend   = max(All_lat);
lonstart = min(All_lon);
lonend   = max(All_lon);

    if ispc 
             slashx = '\';
        elseif isunix 
            slashx = '/';
    end 
   
   
    
    Year_fld = num2str(year); 
    
    if doy < 10 

        DOY_folder = ['00',num2str(doy)];

    elseif   doy >=10 && doy <100

        DOY_folder = ['0',num2str(doy)];

    elseif doy >= 100

         DOY_folder = num2str(doy);

    end 
    
    
     % what type of OS are we running on
    Common_str = ['*',DOY_folder,'-',Year_fld,'*'];
    
    
    %find all the CMN data files to be converted to binary
    datafiles = dir([Data_location,slashx,Common_str]);
    
   
    Totalfiles = length(datafiles);
 
  
    openedfiles = 1;
    numbstatni = 0;
   
    inuse= -1;
    
    while ((openedfiles <= Totalfiles)&&( numbstatni <= numbstatn))
        

        % Get the file i and try to open it 
        filei = datafiles(openedfiles).name;
        fid1 = fopen ([Data_location,slashx,filei],'r');
        filebin =  [Data_location,slashx,filei(1:4),Year_fld,DOY_folder,'.bin'];
   
       
    if exist( filebin, 'file')~=213
      
       
     
        if  fid1 ~= -1
            
            % get the Rx position altitude above sea level in meters,
            % Lat, Log  in degerees
            % located on the 3rd line 
            m = 1;
            while 1 % while loop to get receiver position 
                line = fgetl(fid1);
            
                if ~ischar(line)
                    Rx = [ ];
                    break
                end
                if m>3, break, end
                if (m==3)
                    [t] = sscanf(line, '%f %f %f');
                    Rx = [t(1) t(2) t(3) ];

                end
                m=m+1;
            end % end while loop to get receiver position 
            fclose(fid1); 

            % if we dont have Reciever positions, do not proceed with this file. 
        
            if  ~ isempty (Rx)
                rlat = Rx(1);
                rlon = Rx(2);
                % Is this station within the specified region ?
                % if not neglect
                 if  ( rlat>= latstart) && (rlat<= latend)&& (rlon >= lonstart) && (rlon <= lonend)

                     
                                if  numbstatni > 0
                                    [inuse]=checkmindistance ( Rx, Rxo,min_dist_stn);
                                end 

                                if inuse==1 || numbstatni ==0 

                                    fid1 = fopen ([Data_location,slashx,filei],'r');

                                    out = textscan (fid1, '%f %f %f %f %f %f %f %f %f %f', 'HeaderLines',5);
                                    fclose(fid1);
                                    % check  if there is more than one data sample to
                                    % write 
                                    if length(out {1}) > 1
                                        fid2 = fopen ([Data_location,slashx,filei(1:4),Year_fld,DOY_folder,'.bin'],'w');

                                        if fid2 ~= -1
                                           % if the bin was successfully open, append it
                                           % with  Reciever info
                                            fwrite(fid2,  Rx ,'*double');

                                            % then also add the required data  
                                            % This data will later get read by
                                            % fast_read_CMN4

                                            %  sb  =    JD    Time    PRN    Azm   elev   Lat    Lon   STEC    VTEC    S4   ];  
                                            %  sb  = [out{1} out{2} out{3} out{4} out{5} out{6} out{7} out{8} out{9} out{10}];  
                                            % Lets only select the required parameters in
                                            % our analysis 
                                            sb  = [out{2} out{3} out{4} out{5} out{8} out{9}];
                                            fwrite(fid2, sb ,'*double');
                                            fclose(fid2);
                            
                                            clear out  sb

                                            disp(['Number of cmn station files coverted to binary: ',num2str(openedfiles), 'out of ', num2str(Totalfiles)])
                                            
                                            
                                            Rxo =Rx; % set the previous station to the new station
                                            numbstatni=numbstatni+1;   % increase the station counter 
                                            inuse = -1; % reset the inuse flag 
                                            
                                        else 
                                            clear out  sb
                                            disp(' Please double check the Datalocation path file; failing to create binary file ')
                                        end 


                                    else 
                                      disp([Data_location,slashx,filei, ' has no data!'])
                                    end 
                                else 

                                end 
                 else
                     fclose all;
                 end
                 
           end % No Rx position 
        else
            disp([Data_location,slashx,filei, ' is broken or not avaliable'])
        end
    else
        fclose all ;
        disp([filebin,' already exist'])
    end 
     openedfiles = openedfiles+1;
    end % while end for big I
    disp('Finished coverting *.CMN data to binary format for faster reading')

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
