function [RXo_frm_EOF] = compute_intial_using_EOFs(year,DOY,Time,time_to_collect_rays,Elevation,numbstatn,satheight,Glon,Glat,Gheight,Data_location,length_data_EOF,min_dist_stn)
    % By Nicholas Ssessanga while at Chungnam National University
    % 2019-july-02
dn = DOY; 
[U,No_EOF,All_lat,All_lon,All_HGT] = get_intial_guess_EOFs(Glon,Glat,Gheight,year,dn,Time,length_data_EOF); 

%[H,STEC,~,~]=creat_H_b_bError_anistro(year,DOY, Time,time_to_collect_rays,Elevation,numbstatn,satheight,All_HGT,All_lat,All_lon,Data_location);
[H,STEC,bError,Rxinfo]= creat_H_b_bError_anistro (year,DOY,Time,time_to_collect_rays,Elevation,numbstatn,satheight,...
                                                                                    All_HGT,All_lat,All_lon,min_dist_stn,Data_location);

 [RXo_frm_EOF] = Intial_guess_No(H,STEC,U,No_EOF);

 end 