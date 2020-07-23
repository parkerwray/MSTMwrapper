


%{
    Find optimal weights for making spectral filter out of particles.
%}


clc
clearvars Escat Escat0 Escat180
ff = 0.3;
%[~,idxff] = (min(abs(fill_fractions - ff)));
idxff = 1;
ldamin = 300;
ldamax = 800;
[~,idxldamin] = min(abs(wavelengths-ldamin));
[~,idxldamax] = min(abs(wavelengths-ldamax));
idxlda = idxldamin:idxldamax;

I = cat(3,Ipar,Iper);
Istats = get_statistics(I,3);

qs = squeeze(qsistats.mean(:,idxff,:,:)+qsdstats.mean(:,idxff,:,:));
qa = squeeze(qastats.mean(:,idxff,:,:));

Eabs = (qa(:,idxlda)).';
Escat = (qs(:,idxlda)).';
Escat0 = (squeeze(Istats.mean(:,idxff,:,idxlda,1))).';
Escat180 = (squeeze(Istats.mean(:,idxff,:,idxlda,2))).';
FBR = Escat0./Escat180;
%%



Rideal = zeros(size(wavelengths));
[idxmin, idxmax] = get_box(300, 400, wavelengths);
Rideal(idxmin:idxmax) = 0.0;
[idxmin, idxmax] = get_box(400, 500, wavelengths);
Rideal(idxmin:idxmax) = 0.25;
[idxmin, idxmax] = get_box(500, 600, wavelengths);
Rideal(idxmin:idxmax) = 0.5;
[idxmin, idxmax] = get_box(600, 700, wavelengths);
Rideal(idxmin:idxmax) = 0.25;
[idxmin, idxmax] = get_box(700, 800, wavelengths);
Rideal(idxmin:idxmax) = 0.0;
% [idxmin, idxmax] = get_box(600, 700, wavelengths);
% Rideal(idxmin:idxmax) = 1;

%R = @(w)norm(get_R(w,ff,Escat, Escat0, Escat180)+get_A(w,ff,Eabs)).^2;
R = @(w)norm(abs(Rideal-get_R(w,ff,Escat, Escat0, Escat180))).^2;
%nonlcon = @(w)PMF_constraints(w);

w0 = ones(size(center_radiis))./length(center_radiis);
lb = zeros(size(w0));
ub = ones(size(w0));
Aeq = ones(size(w0));
beq = 1;

w = fmincon(R,w0,[],[],Aeq,beq,lb,ub);


Ropt = get_R(w,ff,Escat, Escat0, Escat180);
norm(Ropt).^2;

%%

figure, 
subplot(2,3,1)
hold on 
bar(center_radiis, w0)
bar(center_radiis, w)
hold off
legend('Original', 'Optimized')
title('Distribution')
lda = wavelengths(idxlda);


subplot(2,3,2)
hold on 
for idx = 1:length(center_radiis)
    plot(lda, log10(FBR(:,idx)))
end
hold off
xlim([min(lda), max(lda)])
title('FBRs')

subplot(2,3,3)
hold on 
for idx = 1:length(center_radiis)
    plot(lda, Escat(:,idx))
end
hold off
xlim([min(lda), max(lda)])
title('Scattering Efficiency')

subplot(2,3,4)
hold on 
for idx = 1:length(center_radiis)
    plot(lda, Eabs(:,idx))
end
hold off
xlim([min(lda), max(lda)])
title('Absorption Efficiency')


Rorig = get_R(w0,ff,Escat, Escat0, Escat180);
subplot(2,3,5)
plot(lda, 100.*Rorig)
hold on 
plot(lda, 100.*Ropt)
plot(lda, 100.*Rideal)
hold off
legend('Original','Optimized','Ideal')
xlim([min(lda), max(lda)])
title('Reflection')


Aorig = get_A(w0,ff,Eabs);
Aopt = get_A(w,ff,Eabs);
subplot(2,3,6)
plot(lda, 100.*Aorig)
hold on 
plot(lda, 100.*Aopt)
hold off
legend('Original','Optimized')
xlim([min(lda), max(lda)])
title('Aborption')







%%
function [cineq, ceq] = PMF_constraints(w)
    
    c_inequality = [];
    c_equality = sum(w)-1;

end

function A = get_A(w,ff,Eabs)

    A = ff.*(w*(Eabs.'));


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


function [idxmin, idxmax] = get_box(ldamin, ldamax, wavelengths)


[~,idxmin] = min(abs(wavelengths-ldamin));
[~,idxmax] = min(abs(wavelengths-ldamax));

end





































