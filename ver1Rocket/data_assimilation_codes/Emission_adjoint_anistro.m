function [EB_djoint] = Emission_adjoint_anistro(modelm, end_of_release,start_of_release,diffusiontimex,everyreleasex,XK,temp_folder) 

% derived symbolic expression for the Emission adjoint 
% --2.2e-15*(xk1*exp(-r1^2/(4*D*k)))/(4*D*pi*k)^(3/2)


EB_djoint =zeros(size(XK(:)));

everyrelease = everyreleasex;

dummydiffusion = start_of_release;
%setofvalues = 0:everyrelease:end_of_release;
cnx=0;

for tp = start_of_release:everyrelease:end_of_release
 
    if tp > diffusiontimex 
         break
    end   

if dummydiffusion >= start_of_release
        
    savetime = roundn(tp,-4);   

    [~,dist_from_source] = get_diffusion_height_anistro (temp_folder,savetime);      


    % get the diffusion time but subtract off the elapsed every 10 second release
    
    t =(diffusiontimex-(everyrelease)*cnx)*60;

    % generate new diffusion mid point between rocket and n(r,t) point
    % this part was computationaly expensive.. 
    % so I ran and saved the data before hand. 
    % this code extracts the saved data which makes the code run fast

    [D1,D2]=get_diffussion_anistro(temp_folder,savetime);

    % compute the adjoint dX/dEo 
    if strcmp(modelm,'H2O')
        % A1 = 1./((D1.*(4*pi*t)).^1.5);
        % 
        % B1 = exp((-1*dist_from_source)./(D1.*(4*t)));
        % 
        % EB_djoint1 =  (-2.1e-15*((A1.*B1).*XK));

        [EB_djoint1]= subadjoint(D1,t,dist_from_source,XK);

    elseif strcmp(modelm,'H2')
        % A2 = 1./((D2.*(4*pi*t)).^1.5);
        % 
        % B2 = exp((-1*dist_from_source)./(D2.*(4*t)));
        % 
        % EB_djoint1 =  (-2.1e-15*((A2.*B2).*XK));

         [EB_djoint1]= subadjoint(D2,t,dist_from_source,XK);

    end 

    EB_djoint = EB_djoint + EB_djoint1;
    
    cnx=cnx+1;



end 

    dummydiffusion = dummydiffusion + everyrelease;
    
end
EB_djoint = EB_djoint;

function [EB_djointi]= subadjoint(Di,t,r2i,Xi)
% Ai = 1./((Di.*(4*pi*t)).^1.5);
% Bi = exp((-1*r2i)./(Di.*(4*t)));
% EB_djointi =  ((-1*2.1e-15)*((Ai.*Bi).*Xi));
const1 = 4*pi*t; 
bot1 = (Di*const1).^1.5;
bot2 = (Di*(4*t));
Ai = 1./bot1; 
expAi = exp((-1*r2i)./bot2);
EB_djointi = ((-1*2.1e-15)*(Ai.*expAi)).*Xi.*150;


