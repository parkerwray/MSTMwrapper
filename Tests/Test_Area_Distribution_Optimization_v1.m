


%{
    Find optimal weights for making spectral filter out of particles.
%}




ff = 0.3;
Escat = [1,1,1,1;1,1,1,1];
Escat0 = [100,0,0,0;0,0,0,100];
Escat180 = [1,1,1,1;1,1,1,1];


R = @(w)norm(get_R(w,ff,Escat, Escat0, Escat180)).^2;
%nonlcon = @(w)PMF_constraints(w);

w0 = [0.1, 0.2, 0.3, 0.4];
lb = zeros(size(w0));
ub = ones(size(w0));
Aeq = ones(size(w0));
beq = 1;

w = fmincon(R,w0,[],[],Aeq,beq,lb,ub);
figure, 
bar(w)
Ropt = get_R(w,ff,Escat, Escat0, Escat180)
norm(Ropt).^2


function [cineq, ceq] = PMF_constraints(w)
    
    c_inequality = [];
    c_equality = sum(w)-1;

end


function R = get_R(w,ff,Escat, Escat0, Escat180)

    sca = w*(Escat.');
    sca0 = w*(Escat0.');
    sca180 = w*(Escat180.');
    
    R = (ff.*sca)./(1+(sca0./sca180));

end



function R = get_Rnorm(R)
    R = norm(R).^2;
end














































