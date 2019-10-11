function [U,No_EOF,All_lat,All_lon,All_HGT] = get_intial_guess_EOFs(Glon,Glat,Gheight,year,dn,hry,length_data)
dn_basket = zeros(1,(length_data*2)+1);
hry_basket = zeros(1,(length_data*2)+1);
year_basket = zeros(1,(length_data*2)+1);

posi=1;
for N =-length_data:length_data
dn_basket(posi)= (dn-N);
hry_basket(posi) = (hry-N);
year_basket(posi) =(year-N);
posi=posi+1;
end 

dn_basket(dn_basket<0) = dn_basket((dn_basket<0))+365; 
dn_basket(dn_basket>365) = dn_basket(dn_basket>365)-365; 



hry_basket(hry_basket<0) = hry_basket(hry_basket<0)+ 24 ; 
hry_basket(hry_basket>24) = hry_basket(hry_basket>24) - 24 ; 



[combined_elecX,~,All_lat,All_lon,All_HGT]=density (Glon,Glat,Gheight,year_basket,dn_basket,hry_basket);

[U,No_EOF] = find_U(combined_elecX);


end



function [U,No_EOF] = find_U(combined_elecX)

%     lt = All_lat; 
%     ln = All_lon; 
%     ht = All_HGT;
% 
%     M = length(lt)*length(ln)*length(ht);
%     N =length(year)*length(dn)*length(hry);
%     data_no_mean = reshape(combined_elecX,M,N);

    [U,D,~] = svd(combined_elecX,'econ');
    % lets find the number of EoFs with most power above 98.8 % 
   % clear V 
    for dd = 1:length(D)
       if  trace(D(1:dd,1:dd)/trace(D))*100 > 99   
          No_EOF = dd; 
           break 
       end 

    end

    if No_EOF > 6 
        No_EOF = 6;
    end 

    U = U(:,1: No_EOF);
end




