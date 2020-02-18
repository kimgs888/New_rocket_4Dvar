    % This is the main code in estimating the North korean ROcket Emissions 
    %  Details can be found in a paper by Ssessanga et al 2018
    %  The 4D?var Estimation of North Korean Rocket Exhaust Emissions Into
    %  the Ionosphere ( Kwangmyongsong?4  Rocket analysis)
    % By Nicholas Ssessanga while at Space Science Lab, Chungnam National University 


    
    
    % before starting, clear and close all currenly open variables. 
    % This will f0.ree up memory 

    clear all ;
    clc;
    close all ;
    fclose all;




    % define the year and DOY being processed 

    Year = 2016; 
    DOYS = 38;   

    % Define your working directory 
    working_dir = '/root/Desktop/Gwangsu/ROCKET5';
    cd (working_dir)
    addpath( working_dir)
    % Define the Location of Data to be assimilated 
    Data_STEC_CMN_folder = '/root/Desktop/4DVar/Suin/Four_Dvar_Nick/Data_STEC_CMN_folder';

    % Define where the final solution should be saved
    Final_solution_flder = '/root/Desktop/Gwangsu/ROCKET5/Final_solution'; 

 
    %                           TEMP FOLDER 
    % all data corresponding to the different time segments within a single assimilation 
    % window are processed and stored in this temp_folder, for quick access during the iteration
    % process. This process highly reducess the computation time!
    % Note: make sure the location of the temp_folder is on a drive with enough free space ~ 1 GB 


    Temp_fld_Loc = '/root/Desktop/Gwangsu/ROCKET5'; 

    %Define minimum Elevation Angle for GPS rays %%
     Elevation=30;

    % Define maximum number of used stations to be used in analysis %% 
     numbstatn=200;
     % 
     min_dist_stn = 0; 
     
     type_cov ='diagonal';
     
    % covergence  parameters 

     B_dump_factor  = 80;%/10; % .5;%(1e-1);% Dump factor on B 
     B_dump_min=1e-23;
     max_while_iter = 80;    % maximum number of iteration that should be performed while 
                             % while inside a while loop 

     chitol= 0.05;           % chi-squared tolorance 
     chivalue=.5;            % chi-squared  value  

    % define length of assimilation window and correlation time (Tau)

     Len_Assim_win = 40;% 15;   % assimilation window length in minutes 

     tauhours = 1;              % correlation time in hours (Tau)

     tauhoursBx =24;       % time to elapse before recomputing new B matrix from background(IRI): 
                           % currently equal to tau but can be changed to
                           % another value as desired

     tauhoursXb =1;        % time to elapse before recomputing new Xb from background(IRI): 
                           % currently equal to tau but can be changed to
                           % another value as desired

     time_sample = 2.5/2;  % sample or length of segments in assimilation window 
                           %(in minutes,  currently 2.5 min). We divide by two because 
                           % the window is centred at the middle value (\..!../)

     time_to_collect_rays = 2.5/2; % same as time_sample above but can be changed 


     increase_sd = 0;      % factor to multiply the data variance 
                           % if the data is less trustworth; variance(Rij, i=j) = (increase_sd*Rij)

     chivalue_rerun=0;     % flag to check whether we should rerun to obtain a better chi_value
     
     
    
    %                   SET UP EMISIONS 
%     
%     a = [1e23 2e23 3e23 3.5e23 4e23 5e23 ];
%     b = 5e23;
%     lenga = length(a);
%     for k = 1 : lenga  
%         QH20  = a(k).^2; 
%         QH2  = b.^2;
    EbH20 =4e23;            % H2O Emissions value used in estimating variance 
    EbH2  =4e23;            % H2 Emissions value used in estimating variance

    QH20  = EbH20^2;       % Variance estimate of the H2O Emisions 
     QH2  = EbH2^2;        % Variance estimate of the H2  Emisions 
    
    EoH20 =0e27;            % Intial background estimated of Hydrogen (H2O) 
    EoH2 = 0e27;            % Intial background estimated of Hydrogen (H2) 
    
    Eb_dif = [];            % Set to store estimated emissions to covergence 
                            % store the intial estimates 
    Eb_dif = [Eb_dif; [( EoH20+EoH2 ) EoH2  EoH20]];
    
    sum_EB_adjnt_X_lambdaH20=0; % Intialize sum of adjoints for H20
    sum_EB_adjnt_X_lambdaH2=0;  % Intialize sum of adjoints for H2
    
    factorEb =10;                % Dump factor for the Emissions 
    Bmatfactor=3.44e-6; 
    
    
    %     Define grid in Quest: 

    % Height 
    Gheight{1}=100:20:450;
    Gheight{2}= 600:200:1400;

    % Latitude 
    latbegin =24;    % Start of Latitude 
    latend =44;      % End Latitude 
    latstep =1;      % Lat resolution  

    % Longitude 
    longbegin=109; % Start of Longitude
    longend= 130;    % End Longitude
    longstep =1;     % Long resolution   

    % Rocket bearing  from the  Sohae space center (124.7?E, 39.6?N)
    Long1 =124.7;
    Lat1  = 40;
    bearig = 180;
    
    min_iter = 5;

    % Assumed rate of Emission effusion along the trajectory
    % we have set this to 10 seconds 
    everyrelease = 10/60;   % In minutes 
    
    % Start and end of Emission release in minutes  
    % These values were inherited from analyses of  previous launches (Taepodong?2). 
    % From the works of Ozeki and Heki (2010), the first engine jettison 
    % completed below 100 km in altitude, at about 2.6 min after launch.
    % Therefore, focus here  will be on the second stage which traversed the ionosphere to a height of 385 km.
    % see The 4D?var Estimation of North Korean Rocket Exhaust Emissions
    % Into the Ionosphere paper
    end_of_release = 6.5; % In minutes 
    start_of_release=2.5; % In minutes 
    
    % Total time of effusion 
    Total_release_time =  (end_of_release-start_of_release);  % In minutes 
    
    
    % parameter used in computing the Background while using EOFs 
    length_data_EOF = 3;
    chibreakicrease=0;
   
%      RMSEK.STRT = 0;
     RMSEK1.STRT = 0;
    % ~~~~~~~~~~~~~~~~~~~~~~  SHOULD NOT EDIT BELOW UNLESS IF KNOW WHAT YOU ARE DOING ~~~~~~~~~~~~~~~~~~~~~~~~

    % Find out whether the platform is pc or unix 
    if ispc 
         slashx = '\';
    elseif isunix 
         slashx = '/';
    end 

    % Add these paths where the routines to be used in computation are located 
    
    %location of the geometry code 
    addpath([working_dir,slashx,'geometry_bath']); 
    %location of of assimilation codes 
    addpath([working_dir,slashx,'data_assimilation_codes']);
    % location of IGRF % goemagnetic code
    % Note: this should be revised and replaced with the corrected version
    addpath ([working_dir,slashx,'IGRF2'])
   % Create_raytrace_mex(working_dir) 
    
    % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                  spacefxn 

    
    
    
    
    
    
    %  ~~~~~~~  set up the temp_folder  ~~~~~~~  
    disp('** 4D-var temp folder setup **')
    locationx =[Temp_fld_Loc,slashx,'Fdvar_temp_folder'];


    % Check and see if this folder does exist. 
    % if this is the first time running at this location create the temp folder

    if exist (locationx, 'dir')
        disp('4D-var_temp_folder already existing in ....')
        disp(Temp_fld_Loc)
        disp ('Only a clean up will be performed!')

        % perform the clean up now 
        tic 
        disp ('cleaning up .......')
        filex = dir([locationx,slashx ,'*.mat']); 
        for f=1:length(filex)
        delete([locationx,slashx ,filex(f).name])
        end 
        disp ('finished cleaning up .......')
        toc 
        
    else % this folder does not exist, lets create a new one 
        [s,mess,messid] = mkdir (Temp_fld_Loc,'Fdvar_temp_folder'); 

        if s~=1 % status for scussefully creating folder, 1 = okay otherwise throw an error
            disp('failed to create Fdvar_temp_folder!')
            error('check temp folder path ..... ')
        end 
    end

    
    % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                  spacefxn 
                  
    %    ~~~~~~~  B matrix folder set up   ~~~~~~~  
    %  define where the B matrix  should be saved every time we generate a new one 
    
    disp('** B matrix temp folder setup **')
    Blocation = [Temp_fld_Loc,slashx ,'Bmatrix_flder'];

    if exist (Blocation, 'dir')
        disp('B matrix folder already existing in ....')
        disp(Temp_fld_Loc)
        disp ('Only a clean up will be performed!')
        % perform the clean up now of P matrix folder
        tic 
        disp ('cleaning up B matrix folder.......')
        filex = dir([Blocation,slashx,'*.mat']); 
        for f=1:length(filex)
        delete([Blocation,slashx,filex(f).name])
        end 
        disp ('finished cleaning up  B matrix folder .......')
        toc 

    else
        [s,mess,messid] = mkdir (Temp_fld_Loc,'Bmatrix_flder'); 
        if s == 1 
            disp('succefully created B matrix folder....')

        else % status for scussefully creating folder, 1 = okay otherwise throw an error
            disp('failed to create B matrix folder!')
            error('check temp folder path ..... ')
        end 
    end

   % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                  spacefxn 
                  
                  
   % ~~~~~~~  Final solution folder set up   ~~~~~~~ 

    disp('** Final solution folder set up **')

    if exist (Final_solution_flder, 'dir')
        disp('Final solution folder already existing in ....')
        disp(Final_solution_flder)    
    else
        [s,mess,messid] = mkdir (Final_solution_flder); 
        if s == 1 
            disp('succefully created Final solution folder....')

        else % status for scussefully creating folder, 1 = okay otherwise throw an error
            disp('failed to create Final solution folder!')
            error('check Final solution folder path ..... ')
        end 
    end
    
    
     %  ~~~~~~~       end of SET UP  ~~~~~~~  % 
 
     

    % Grid specifications 
    disp ('** Grid specifications **')

    Grdlatst =latbegin; 
    Grdlatend =latend;
    Grdlongst =longbegin;
    Grdlongend = longend;

    Glat{1}= latbegin:latstep:latend;
    Glon{1}=longbegin:longstep:longend;

    [All_lat,All_lon,All_HGT]=gridvectors(Glon,Glat,Gheight);



    disp(['Number of Latitudes: ... ',num2str(length(All_lat(:))) ])
    disp(['Number of Longitudes: ... ',num2str(length(All_lon(:))) ])
    disp(['Number of Heights: ... ', num2str(length(All_HGT(:)))])
    numgrid = length(All_lat(:))* length(All_lon(:))* length(All_HGT(:));
    disp(['Total Grid cells: ', num2str(numgrid)])

    % assume satellites are located at this altitude 
    satheight = All_HGT(end)*1E3+50E3;
    % Grid height boundary
    hgtbgn = min(All_HGT);
    hgtend= max(All_HGT);
    
    % Define where the STEC data to be ingested is located   
    % This data should be in .CMN format, processed using the Gopi program 
      
  
     % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                   spacefxn 

     % generate the mesh to be used in computing defusion at different
     %  heights
     
     disp('Generating Mesh and trajectory to be used in computing Difusion')
     [xi, yi, zi] = meshgrid(All_lat,All_HGT,All_lon);
     % Call this function to generate the rocket Trajectory assumed from
     % the work of  Ozeki and Heki (2010), but with a southward direction. 

     [Lat2,Long2,Heighti, Traje_time]=Rocket_Trajectory_frm_Hekie_paper(Long1,Lat1,bearig);
     disp('Finished generating Mesh and Trajectory')
     
     disp('Generating and saving Difusion for different time steps before computation')
     % To approximate a more realistic anisotropic diffusion that depends on altitude, 
     % model 1 presented in Furuya and Heki (2008) was adopted:
     % Computation cost is highly reduced is this diffusion is computed outside the estimation loop 
     
     save_difused_nrt_mendillo_anisttrotopic_new(xi, yi, zi,Traje_time,Lat2,Long2,Heighti,hgtbgn,hgtend, Total_release_time,...
                                                                 everyrelease, end_of_release, start_of_release, locationx)
     disp('Finished Generating and saving Difusion to be used later in Computation')
     
     
     % set to store all chi estimated values 
     chi_set = []; 
  
for DOY =DOYS
    % Rocket was launched 00:30 UT and the second stage which is of
    % interest here started ~ 2.5 minutes after launch 
    
    Time_of_start =0.5 + start_of_release/60;    % Analysis Time 4Dvar should Start the processing 
 

    start_assim_win = Time_of_start;             % set the begining of the assimilation win to Time_of_start 

    % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                   spacefxn 
                   
    % Define where the STEC data to be ingested is located 
    % this data should be in .CMN format, processed using the Gopi program 
    disp('** Data folder check **')
    
   [Data_location]= fomulate_data_location_flder (Data_STEC_CMN_folder,slashx,Year,DOY);
    
    % Covert the available data into binary format 
    % This data format speeds up the file reading procedure 
     covert_cmn_data2binaryx ( Data_location,Year,DOY,All_lat,All_lon,numbstatn,min_dist_stn) 

   % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                      spacefxn 
    

     % create the assimilation window container 
     % Last time sample will be start point + length of window 
     % remember this time is in hours so covert every thing to hours

     End_of_Assim_win =  start_assim_win + Len_Assim_win/60;%hours 
 

     %  Vxi vessel contains all  the time segments within the assimilation
     %  window moving backwards in time. This will be the first assimilation
     %  window to be processed. The next window will be updated at the end of
     %  the first while loop. 
     Vxi = End_of_Assim_win :-2*(time_sample/60):(start_assim_win);


    %Intialize flags to be used in computation 
     LastcounterXs=1;
    % check whether to compute new Background 
     corrletionT = 0; 
     corrletionTx =0;
    % check whether to compute new B matrix
     corrBcompuT =0;
     corrBcomputx = 0;
    % chi break flag 
     chibreakicrease=0;
  

    % Set this value which will help to reset the dump factor 
    % if its modified with the code. For every new assimilation window, 
    % this factor should be restored to the original value.

    B_dump_factoro = B_dump_factor;

    
    % intialize structure to store data 
    xcount=1;
    Xanalysis.start = 0;
    X.start = 0;
    cold_run =1;
    usingfiledata=1;
    simulation_analysis = 0;
    chidiffset = [];
    

    %  Setup is done. Lets begin the computation 
    while 1    % this while loop goes through all the Time samples to the End
               % this is the bigger Loop start time ----> End_time


        

        sum_EB_adjnt_X_lambdaH2=0;
        sum_EB_adjnt_X_lambdaH20=0;
        iterations = 1;
        gettingData = 1;
        countwhile =1;
        count_B_dump_flag = 1;
        count_B_dump=1;

        % This is the Assimilation Window optmisation 
        while 1
      
            % we must work moving backwards in time because of the adjoint operation
            
            intial = 1;
            timecounterXs = LastcounterXs; 
            timecounter = 1;
            chi_better_b4_iter = 0; 

            % get the number samples in the assimilation window 
            number_of_time_samples = length(Vxi); 

            % call this fxn to create some space 
            spacefxn 

            disp ('** Time assimilation details **')
            disp (['Length of assimilation Window,...', num2str(Len_Assim_win), ' min'])
            disp (['Number of sumples in assimilation window,...', num2str(number_of_time_samples)])
            disp (['Length of each sample,...', num2str(2*time_sample), ' min'])

            % setup just too store number of GPS data in assimilation 
            GPS_data_count =[];
            clear STECRxPRN
            STECRxPRN.srt = 0;


            % Lets go through the time samples while moving back in time


            for simulationtime =Vxi    
               

                % this is the current time to analyse along the time trajectory 
                Time =simulationtime;


                % processing time in minutes
                Total_process_time  =  (Time - start_assim_win)*60;
             
               


                %******************** data ****************************************% 
                % This only to check if we using data from file(obervation ) or simulation 
                % but currently this is set to using observed data 

                if usingfiledata==1 


                % if this is the first iteration collect all data and store it
                % Here we store the STEC, H, Rx and Tx 
                % We will not have to look for this data again once collected.

                    if iterations ==1 && (gettingData == 1)&&(countwhile==1)&&(chivalue_rerun==0)


                        Current_tm_strhr = floor(Time);
                        Current_tm_strmin = 60* (Time-Current_tm_strhr);

                        if Current_tm_strhr<10 
                            Current_tm_strhrs = ['0' num2str(Current_tm_strhr)];  
                        else 
                              Current_tm_strhrs = num2str(Current_tm_strhr); 
                        end 

                        if Current_tm_strmin <10 
                            Current_tm_strmins = ['0' num2str(Current_tm_strmin)];  
                        else 
                            Current_tm_strmins = num2str(Current_tm_strmin); 
                        end 


                        disp (['Epoch being processed ', num2str(Year), '_',num2str(DOY),'_',Current_tm_strhrs,'h',Current_tm_strmins,'min'])
                        

                        % ~~~~~~~  call this fxn to create some space  ~~~~~~~  
                                              spacefxn 
    

                        disp('Getting data in 1st iteration ...... ')
                        
                        % This function gets the H, b = STEC, and corresponding b Error. 
                         [H,STEC,bError,Rxinfo]= creat_H_b_bError_anistro (Year,DOY,Time,time_to_collect_rays,Elevation,numbstatn,satheight,...
                                                                                    All_HGT,All_lat,All_lon,min_dist_stn,Data_location);

                        fclose all;
                   
                     
                     

                    else
                        % if all the data with in the assimilation is stored in
                        % the next assimilation we will not have to search for the data
                        % and create H. we just have to use the data in the temp folder

                        [H,STEC,bError] = loadsavedHSTEC_temp(Total_process_time,simulation_analysis,time_to_collect_rays,time_sample,locationx);
                    %    disp('Fast computing iterations now ...... ')

                    end 

                end 

                %******************** data ****************************************% 






                %********************************* Xb, and B_matrix **************************%
                % Here we compute the background and its corresponding Covariance
                % matrix 
                % only generate these quantities if  its the intial iteration with in
                % the  assimilation window. other wise this value is the same 
                if intial && iterations == 1


                    % compute the Xb densities using the IRI 2016 module  
                    % and store them in RXo
                    if corrletionT ~= 1

%                           if (chivalue_rerun==0)
%                                [~,RXb,All_lat,All_lon,All_HGT]=density (Glon,Glat,Gheight,Year,DOY,Time);
%                                 RXbiri = RXb(:);
%                               length_data_EOF = 2;
%                               time_to_collect_raysx =10;
%                                load('EOFXo.mat')
%                                RXb4 = Xo;
%                              [RXo_frm_EOF] = compute_intial_using_EOFs(Year,DOY,Time,time_to_collect_raysx,Elevation,numbstatn,...
%                                                                   satheight,Glon,Glat,Gheight,Data_location,length_data_EOF,min_dist_stn);                              
%                                                                RXb=RXo_frm_EOF(:);
%                                                               Xo = RXb;
%                                                              
%                  save EOFXo.mat Xo
                                                             load('EOFXo.mat')
                                                             RXb = Xo;
%                                                               plot(Xo,'*-b')
%                                                               hold on
%                                                               grid on
%                                                               plot(RXb4,'o-r')
%                                                               
%                                                               pause
% % %                              
                             
% 
%                           else
%                                RXb = RXb(:); 
%                           end 
                        
                          if cold_run ==1 
                                 RXo=RXb(:);
                                 
                                 
                                 cold_run =0;
                          end 
                       

                        % set the flag to compute Background to 1, to indicate that its
                        % Xb is already generate
                        corrletionT =1;
                    end 
                    % compute the B using the IRI 2016 module  
                    % and store file in Plocation folder
                    if corrBcompuT ~=1
                         % set the flag to compute Background coveriance (B)to 1, to indicate that its
                         % already generated
                       
                           
                         
                    
                       
                       cov_B_matrix_anistro(All_lat,All_lon,All_HGT,RXo,Blocation,Bmatfactor,type_cov) 
                       corrBcompuT=1;
                       
                    end

                     fclose all;

                end 

                %******************** Xb B_matrix ****************************************%



                %******************** Xbi ****************************************%
                % Set intial for all time samples as the background
                % Here we are assuming that within 15 minutes the background will not
                % change that much. so lets set intial Xn to Background. 
            
                if iterations == 1 && intial 
                          
                   
                    [~,~,Xo,X] = set_all_background(X,Vxi,RXo,timecounterXs) ;
                    
                   
                   

               

                end 
                %******************** Xb ****************************************%



                %******************** Forecast ****************************************%

                % propagation_time = Total_process_time*60;% time to propagate the density 
                diffusiontimex  = Total_process_time;
                % [XK]=forecast(Xo,RXo,tau,deltaT,propagation_time);%(Xo,tau,deltaT,propagation_time); 
                [XK,hole ] = difused_nrt_anistro(end_of_release,start_of_release,diffusiontimex,everyrelease,EoH20,EoH2,Xo,locationx,time_sample);
%                  [RMSEK]=fnxRmse(XK,H,STEC,diffusiontimex,iterations, RMSEK);
                  
                %******************** end forecast ****************************************%

%

                % *****************Quality control *********************************************%
                 if iterations ==1 && (gettingData == 1)&&(countwhile==1)
                      % only perform quality control if the 
                      % data are enough for
                      % for statics 
                      load Tottec.mat
%                       if size(STEC,1) >= 100 
%                          [H,STEC,bError] = Quality_control (H,STEC,XK,bError);
%                       
%                       else 
%                           disp('Not enough data for Quality control statistics')
%                       end 


                    % save the data after Quality control
                    % Save the corresponding data to data (H, STEC and bError) temp folder in accordance to time                                 
                    saveHSTEC_temp(Total_process_time,H,STEC,bError,simulation_analysis,time_to_collect_rays,time_sample,locationx)

                    % Save the corresponding data to data (H, Rx and Tx ) temp folder in accordance to time 
                    saveRxLat_RxLon_PRN_temp(Total_process_time,Rxinfo.RxLat,Rxinfo.RxLon,Rxinfo.prn,Rxinfo.time,...
                                                      simulation_analysis,time_to_collect_rays,time_sample,locationx)
                 end 
                 
                 [RxLat ,RxLon,Time,PRN] = loadsavedRXTX_temp(Total_process_time,simulation_analysis,time_to_collect_rays,time_sample,locationx);
                 [STECRxPRN] = saveRxLat_RxLon_PRN_TEC(Total_process_time,RxLat ,RxLon,PRN,Time, H, XK, STEC,STECRxPRN); 
                 [RMSEK1] = fnxRmse2(XK,H,STEC,diffusiontimex,iterations, RMSEK1,RxLat,RxLon,PRN);
                % *****************Quality control *********************************************%

  




                %************************ Compute lambda_N *************************************%

                if intial   

                        [lambda_k_1,~] = Compute_lambda_Nx(H,STEC,XK,bError,Time,increase_sd );

%                         modelm='H2O';
% 
%                         [sum_EB_adjnt_X_lambda] = sumofadjoints_anistro(modelm,end_of_release,start_of_release,...
%                                                         sum_EB_adjnt_X_lambdaH20,diffusiontimex,lambda_k_1,everyrelease,XK,locationx,time_sample);
% 
%                         sum_EB_adjnt_X_lambdaH20 = sum_EB_adjnt_X_lambda;
% 
%                         modelm='H2';
% 
%                         [sum_EB_adjnt_X_lambda] = sumofadjoints_anistro(modelm,end_of_release,start_of_release,...
%                                                                     sum_EB_adjnt_X_lambdaH2,diffusiontimex,lambda_k_1,everyrelease,XK,locationx,time_sample);
% 
%                         sum_EB_adjnt_X_lambdaH2 = sum_EB_adjnt_X_lambda;


                        intial = 0; 
       
                end 

                %************************end of Compute lambda_N *************************************%




       
                %use lambda_N to compute the next  adjoint (lambda_k) as you move back in time 
             

                if timecounter > 1 && timecounter < length(Vxi)        
                       
                        [lambda_k,~,~] = Compute_lambda_k_anistro(H,STEC,XK,bError,Time,increase_sd,lambda_k_1,end_of_release,...
                                                                                start_of_release,diffusiontimex,everyrelease,EoH20,EoH2,locationx,time_sample);
                                                                           


                        modelm='H2O';
                     
                        [sum_EB_adjnt_X_lambda] = sumofadjoints_anistro(modelm,end_of_release,start_of_release,...
                                                        sum_EB_adjnt_X_lambdaH20,diffusiontimex,lambda_k_1,everyrelease,XK,locationx,time_sample);
                                                  
                        sum_EB_adjnt_X_lambdaH20 = sum_EB_adjnt_X_lambda;

                        modelm='H2';
                       
                        [sum_EB_adjnt_X_lambda] = sumofadjoints_anistro(modelm,end_of_release,start_of_release,...
                                                                    sum_EB_adjnt_X_lambdaH2,diffusiontimex,lambda_k_1,everyrelease,XK,locationx,time_sample);
                                                               



                        sum_EB_adjnt_X_lambdaH2 = sum_EB_adjnt_X_lambda;
                        
                        

                        lambda_k_1= lambda_k;
                        

                end

                
                  if simulationtime == Vxi(end) && iterations ==1 
                       
                      
                          [sd] = R_covariance_Matrix2x(STEC,bError,Time,increase_sd);

                     
                          chivaluenew = sum((((H*Xo-STEC)).^2)./(diag(sd)))/length(STEC);
                         
                       
                          if chivaluenew <=chivalue

                          chibreakicrease=1; 
                          chi_better_b4_iter = 1; 
                          break
                          end
                   end 
                
                
               
                %*******************************************************************

                timecounter = timecounter +1;
                timecounterXs = timecounterXs+1;
                
            end % Assimilation window end 




            % Lets Compute the Xo using Xo = Xb + B*?o; Xb here is previous Xo
            if chi_better_b4_iter ~= 1 
                      % call this fxn to create some space 
                                spacefxn 
    %                 disp ('** computing New solution **')
    %                 disp (['For epoch: ', num2str(Year), '_',num2str(DOY),'_',Current_tm_strhrs,'h',Current_tm_strmins,'min'])
                    prev_Xo = Xo; % store previous Xo before any new adjustments
                                  % if anything went wrong in the computation, then we willset 
                                  % the densities to back to the prevous value 
                    EoH2_prev = EoH2;
                    EoH20_prev = EoH20;
                 
                 
                           % ~~~~~   Perform the computation     ~~~~~~~~ 
                           % Load Bmatrix stored in Blocation
                    B = load([Blocation,slashx,'Bmatrix.mat']);
                    B = (B.B);
                           %            New solution            %
                    if    B_dump_factor < B_dump_min %(abs(diffChi)> 1.5 *chitol && abs(diffChi)) < 3*chitol %||
                          B_dump_factor = B_dump_min;
                     end
                    Xo = Xo + B_dump_factor *( B*lambda_k(1:length(B(1,:))));
               
           
                    % Set a condition to limit any densities from 
                    % going below the minimum value here this value is set as 1.0e6 
                    % Nontheless, analysis densities never go below this value from my exprience 
                    [Xo] = set_min_density(Xo);
                    
                    %update the Emissions 
                    EoH20 =  (EoH20 + factorEb*QH20*sum_EB_adjnt_X_lambdaH20);
 
                    EoH2 = (EoH2 + factorEb*QH2*sum_EB_adjnt_X_lambdaH2);
          
                   
                    
              
                    clear B
                    %           Done with New solution 

                    %              Chi squared value  
                                    spacefxn 
                    disp ('** Determing goodness of fit using Chi test **')
                    % This will help us determine how good the new solution is. 
                    % or The goodness of our model fit and choice of parameters 
                    % See underst_chi.pdf in document folder 
                    % ?_(Xi-?)? (var), here we take ? = mean as the measured STEC. 

                    % so lets get the variance (var)
                
                    [sd]=R_covariance_Matrix2x(STEC,bError,Time,increase_sd);
                  

                    % now compute the chi squared value using the previous estimate and the 
                    % current estimate
                    % The lowe chi sqaured value is the better the solution 
                    chivaluenew = sum((((H*Xo-STEC)).^2)./(diag(sd)))/length(STEC);
                    chivalueold = sum((((H*prev_Xo-STEC)).^2)./(diag(sd)))/length(STEC);

                    diffChi = (chivalueold-chivaluenew);
                    % Store these values to be plotted later. 
                    % remember we assumed the emissions to be deposited
                    % along the trajectory every  10 seconds 
                    % to get the amount per second lest divide by 10 
                    Eb_dif = [Eb_dif; [( EoH20+EoH2 )/10 EoH2/10  EoH20/10] ]
                    disp('****************************************************')
                    Eb_dif2 = diff(Eb_dif)
                   
                    chi_set = [chi_set; chivaluenew];
                    disp(['Old chi value...', num2str(chivalueold)])  
                    disp(['New chi value...', num2str(chivaluenew)]) 
                    disp(['updated H2 Emissions (EH2)...', num2str(EoH2/1e27), 'E26 per sec'])
                    disp(['updated H2O Emissions (EH2O) ...', num2str(EoH20/1e27),'E26 per sec'])
                
                    
                    %  Check for covergence, since the emissions are
                    %  starting from 0 as the initial value, they can only
                    %  increase +ve. 
                    %  if the gradient starts becoming negative, 
                    % terminate: We have reached the optimal value
                    % Plot the Result
                    [meanval] =  Compute_meanval(RMSEK1)
                    disp (RMSEK1.RMSE_7_5)
                    num_iter = iterations
%                     length(chi_set)> min_iter && (meanval(end)-meanval(end-1))>
                     if Eb_dif(end,1) > 5.5e+26 && diffChi < 0 %length(chi_set)> min_iter && (Eb_dif(end,1)-Eb_dif(end-1,1))<0 || chivalueold>20%(meanval(end)-meanval(end-1))>0 %(Eb_dif(end,1)-Eb_dif(end-1,1))<0  && length(chi_set)> min_iter || ((abs(diffChi) <=1e-2)&& length(chi_set)> min_iter)||...
%                                 ((diffChi<1e-2)&& length(chi_set)> min_iter) || iterations >300
                                figure
                                [AX,H1,H2] = plotyy(1:length(Eb_dif(:,1)),Eb_dif,1:length(chi_set),chi_set);
                                set(get(AX(1),'Ylabel'),'String','Emissions per sec') 
                                set(get(AX(2),'Ylabel'),'String','chi sqaured at Xo covergence') 
                                set(get(AX(2),'Xlabel'),'String','Iterations') 
                                set(H1,'LineStyle','--', 'linewidth',2.5)
                                set(H2,'LineStyle',':','linewidth',2.5)
                                legend(H1,'H{_2}O + H_2','H{_2}','H{_2}O')
                                grid on 
                                
%                                  plotrmse
                                 plotrmse1
                                  plotstec_timetest2 (STECRxPRN,Tottec)
                             chibreakicrease=1;
                         break 
                     end 



            else 

                 % call this fxn to create some space 
                                spacefxn 
                 disp('chi better before iterations....')
                 disp(['Chi value...', num2str(chivaluenew)])
                 chivalueold=chivaluenew;
                 chi_better_b4_iter=0;



            end 



            % lets reverse some flags if we are done with iteration
            gettingData=0;
            intial = 1;


            iterations = iterations+1;
            
            if iterations > 800
                 chibreakicrease=1;
                break
            end


        end  

        fclose all;
       
   if chibreakicrease==1;
       break
   end
    end
    cold_run =1;
end
% RMSEK1 = rmfield(RMSEK1,"RMSE_0");
% RMSEK1 = rmfield(RMSEK1,"RMSE_10");
% RMSEK1 = rmfield(RMSEK1,"RMSE_12_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_15");
% RMSEK1 = rmfield(RMSEK1,"RMSE_17_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_20");
% RMSEK1 = rmfield(RMSEK1,"RMSE_22_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_25");
% RMSEK1 = rmfield(RMSEK1,"RMSE_27_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_2_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_7_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_30");
% RMSEK1 = rmfield(RMSEK1,"RMSE_32_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_35");
% RMSEK1 = rmfield(RMSEK1,"RMSE_37_5");
% RMSEK1 = rmfield(RMSEK1,"RMSE_40");
% RMSEK1 = rmfield(RMSEK1,"RMSE_5");
% STECRxPRN = rmfield(STECRxPRN,'TEC_34_47_126_32_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_34_47_126_32_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_34_68_126_9_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_34_68_126_9_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_35_18_126_91_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_35_18_126_91_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_35_28_126_52_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_35_28_126_52_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_35_62_126_97_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_35_62_126_97_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_01_126_76_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_01_126_76_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_46_126_8_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_46_126_8_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_78_126_49_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_78_126_49_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_89_126_82_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_36_89_126_82_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_01_127_27_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_01_127_27_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_22_127_23_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_22_127_23_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_28_127_05_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_28_127_05_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_42_126_69_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_42_126_69_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_75_126_74_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_75_126_74_20');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_9_127_06_15');
% STECRxPRN = rmfield(STECRxPRN,'TEC_37_9_127_06_20');
% clear RMSEK1;clear STECRxPRN;
%     end
