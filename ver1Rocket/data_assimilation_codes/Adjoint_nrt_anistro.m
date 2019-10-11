 function [Adjnt] = Adjoint_nrt_anistro(end_of_release,start_of_release,diffusiontimex,everyreleasex,EoH2O,EoH2,XK,temp_folder)
 % This fuctions computes the Jacobian 

% Evaluate the adjoint
% evaluate  = 
%             d[Xk+1]
% Jacobian  = -------         
%             d[Xk]

everyrelease = everyreleasex;

%diffusiontime =(end_of_release-start_of_release);

nrt =zeros(size(XK(:)));
nrtH2=zeros(size(XK(:)));
cnx=0;

dummydiffusion = start_of_release;
%setofvalues = start_of_release:everyrelease:end_of_release;


for tp = start_of_release:everyrelease:end_of_release

    if tp > diffusiontimex 
         break
    end   



    if dummydiffusion >= start_of_release
        
        savetime = roundn(tp,-4);     
        % [~,indxh] = min(abs(Traje_time-tp)); 

        % extract the the vector containing r^2 ==> dist_from_source
        % note: dist_from_source is sqaured distance if you need real distance
        % please sqrt 

        % [~, dist_from_source ] = new_diffusion_height(Heighti,Lat2,Long2,yi,xi,zi,indxh ) ;               


        [~,dist_from_source] = get_diffusion_height_anistro (temp_folder,savetime);  


        % get the diffusion time but subtract off the elapsed every 10 second release
        t = (diffusiontimex-(everyrelease)*cnx)*60;
        
        % generate new diffusion mid point between rocket and n(r,t) point
        % this part was computationaly expensive.. 
        % so I ran and save the data before hand. 
        % this code extracts the saved data which makes the code run fast
        [D1,D2]=get_diffussion_anistro(temp_folder,savetime);



        %%%%%%%%%%%%%%%% HO2 %%%%%%%%%%%%%%%%
        % A1 = SoHO2./((D1.*(4*pi*t)).^1.5);
        % expA1 = exp((-1*dist_from_source)./(D1.*(4*t)));
        % C1 = (A1.*expA1);

        [C1] = subdifuse (EoH2O,t,D1,dist_from_source);

        %%%%%%%%%%%%%%%%H2%%%%%%%%%%%%%%%%
        % molcule = 'H2'; 

        % A2 = SoH2./(D2.*(4*pi*t)).^1.5;
        % expA2 = exp((-1*dist_from_source)./(D2.*(4*t)));
        % C2 = (A2.*expA2);

        [C2] = subdifuse (EoH2,t,D2,dist_from_source);

        cnx=cnx+1;

        % add up the release effect 
        nrt = nrt + C1;
        nrtH2 = nrtH2 + C2;
      
     
    end 

    dummydiffusion = dummydiffusion + everyrelease;
end


%             d[Xk+1]
% Jacobian  = -------         
%             d[Xk]
    Px = 1 - (((nrt+nrtH2).*(2.2e-15 + 2.0e-15)/2)).*150;

    clear nrt nrtH2 C1 C2

    Adjnt = spdiags(Px ,0,length(Px),length(Px ));

 


function [Ci] = subdifuse (Soi,t,Di,r2i)
% Ai = Soi./((Di.*(4*pi*t)).^1.5);
% expAi = exp((-1*r2i)./(Di.*(4*t)));
% Ci = (Ai.*expAi);
% const1 = 4*pi*t; 
% bot1 = (Di*const1).^1.5;
% bot2 = (Di*(4*t));
% Ai = Soi./bot1; 
% expAi = exp((-1*r2i)./bot2);
% Ci= (Ai.*expAi);
const1 = 4*pi*t; 
bot1 = (Di*const1).^1.5;
bot2 = (Di*(4*t));
Ai = Soi./bot1; 
expAi = exp((-1*r2i)./bot2);
Ci= (Ai.*expAi);




