


%{
    Find optimal weights for making spectral filter out of particles.
%}


clc
clearvars Escat Escat0 Escat180
ff = 0.4;
[~,idxff] = (min(abs(fill_fractions - ff)));

ldamin = 600;
ldamax = 750;
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
%R = @(w)norm(get_R(w,ff,Escat, Escat0, Escat180)+get_A(w,Eabs)).^2;

R = @(w)norm(abs(0.3-get_R(w,ff,Escat, Escat0, Escat180))).^2;
%nonlcon = @(w)PMF_constraints(w);

w0 = [1, 1, 1, 1, 1]./5;
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
    plot(lda, FBR(:,idx))
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
hold off
legend('Original','Optimized')
xlim([min(lda), max(lda)])
title('Reflection')

%%
function [cineq, ceq] = PMF_constraints(w)
    
    c_inequality = [];
    c_equality = sum(w)-1;

end

function A = get_A(w,Eabs)

    A = w*(Eabs.');


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














































