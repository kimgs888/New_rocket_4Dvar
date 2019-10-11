function [sum_EB_adjnt_X_lambda] = sumofadjoints_anistro(modelm,end_of_release,start_of_release,sum_EB_adjnt_X_lambda,diffusiontimex,lambda_k_1,everyreleasex,XK,temp_folder)
                                                                       


if strcmp(modelm,'H2O')
   
[EB_djoint] = Emission_adjoint_anistro(modelm,end_of_release,start_of_release,diffusiontimex,everyreleasex,XK,temp_folder);                                                                       
                                                                     
sum_EB_adjnt_X_lambda = sum_EB_adjnt_X_lambda + EB_djoint'*lambda_k_1;

elseif strcmp(modelm,'H2')
    
[EB_djoint] = Emission_adjoint_anistro(modelm,end_of_release,start_of_release,diffusiontimex,everyreleasex,XK,temp_folder);                                                                       
                                                                     
sum_EB_adjnt_X_lambda = sum_EB_adjnt_X_lambda + EB_djoint'*lambda_k_1;

end 