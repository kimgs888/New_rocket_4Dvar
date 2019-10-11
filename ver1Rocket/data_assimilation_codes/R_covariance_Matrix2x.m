function [R]=R_covariance_Matrix2x(STEC,bError,Time,increase_sd)
    % This function computes the data variance 
    % and stores it in a diagonal matrix R 
    % By Nicholas Ssessanga while at Chungnam National University
    % see the discription in a paper imaging south African ionosphere 
    % using 4Dvar by Nicholas Ssessanga
    % 2018-july-04

    R = sparse(length(STEC),length(STEC));

    for II = 1:length(STEC)
        sigmax = sqrt((STEC(II)*(10/100))^2 + (bError(II))^2) + increase_sd*sqrt((STEC(II)*(10/100))^2 + (bError(II))^2);

        if sigmax > 7e17
          sigmax =7e17;

        end 

        R(II,II)=sigmax^2;

    end

end 

