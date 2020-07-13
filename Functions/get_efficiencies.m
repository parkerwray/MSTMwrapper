function [qe, qa, qsi, qsd] = get_efficiencies(data)

    data = data{1};

    % Grab efficiencties associated with center particle. 
    qe = (data(1).Qext_par+data(1).Qext_per)/2;
    qa = (data(1).Qabs_par+data(1).Qabs_per)/2;
    qsi = (data(1).Qsca_par+data(1).Qsca_per)/2;
    qsd = qe - qa - qsi;
    
end
    
    
    
    






