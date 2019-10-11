function [] = cov_B_matrix_anistro(All_lat,All_lon,All_HGT,RXo,Blocation,Bmatfactor ,type_cov)
% This definition of the Covariance matrix is based on the 
% paper by Bust et la 2004 
% or See Assimilation of Multiple Data types to a Regional Ionosphere Model with a 3D-Var algorithm (IDA4D)
% Chalachew et la., (2019) 

disp('Computing the Covariance Matix ..........')
if strcmp(type_cov,'full')
% calculate the total number of points in defined grid
np = length(All_HGT)*length(All_lat)*length(All_lon);
ras_ind = 1:np; % these will be the indicies of a rasta_scaned vector()
epsilon = 2d-1; % this will be used in defining the cut off point, beyond which every thing is set to zero
				% Zeros make the matrix sparse and hence easy to deal with in computation
  pCorrVal = zeros(np,np);               
				
		%  ~~~~~~~~~~~~ Horrizontal and Vertical Correlations ~~~~~~~~~~~~~~~~~

% height used in computing the horizontal correlation 
All_HGT300 = 300;

% call this fxn to compute the horrizontal correlation 
[hcorrPtr]= ComputeHorizontalCorr (All_lat,All_HGT300,All_lon);
% call this fxn to compute the vertical correlation 
[vcorrPtr ]= ComputeVerticalCorr(All_lat,All_HGT,All_lon);

% reshape the input density into Vertical and Horrizontal shape
% these will be used in computing the assumed Variance values 
sigma_ptr = reshape(RXo(:),length(vcorrPtr),length(hcorrPtr));


% ~~~~~~~~~~~~~~ Generate the Covariance matrix ~~~~~~~~~~~~~~ 
k =1;
  for i = 1:length(hcorrPtr)
      for j = 1:length(vcorrPtr)
         row =vcorrPtr(j,:)'*hcorrPtr(:,i)';
         sigma = sigma_ptr(:,i)*sigma_ptr(j,:);
         sigma = 0.4*sigma(:);
         row = row(:); 
       
         pCorrVal(k, ras_ind(row> epsilon)) =  sparse(row(ras_ind(row> epsilon)).*sigma(row> epsilon));
         
    
   k =k+1;
      end
  end
B = sparse(pCorrVal);
cd(Blocation)
save Bmatrix B 
clear B


elseif strcmp(type_cov,'diagonal')
    
    
    
B = spdiags((RXo.^2).*Bmatfactor ,0,length(RXo ),length(RXo));
cd(Blocation)
save Bmatrix B 
clear B


end 
disp('Finished computing the Covariance Matix ..........')

