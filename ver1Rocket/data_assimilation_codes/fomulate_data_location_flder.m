function [Data_location]= fomulate_data_location_flder (Data_STEC_CMN_folder,slashx,Year,DOY)
 
% This function fomulates the datalocation folder 
% By Nicholas while at CNU 
    Year_fld = num2str(Year); 

    if DOY < 10 

        DOY_folder = ['00',num2str(DOY)];

    elseif   DOY >=10 && DOY <100

        DOY_folder = ['0',num2str(DOY)];

    elseif DOY >= 100

         DOY_folder = num2str(DOY);

    end 

    Data_location = [Data_STEC_CMN_folder,slashx,Year_fld,slashx, DOY_folder];
    
    
     if exist (Data_location, 'dir')
        
     
        numb_filex = length(dir ([Data_location,slashx,'*.cmn']));
        
        if numb_filex ==0
        numb_filex = length(dir ([Data_location,slashx,'*.Cmn']));
        end 


        if  numb_filex > 1
            numb_files = length(dir([Data_location,slashx,'*',DOY_folder,'-',Year_fld,'*']));
            disp ( [num2str(numb_files),  ' CMN data files available for processing'])
        else 
           error('No cmn files available in folder, recheck folder..... ')   
        end 

    else
         disp(Data_location)
         error('data folder does not exist, recheck folder..... ')  
    end