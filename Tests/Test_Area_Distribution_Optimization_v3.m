


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
Rideal(idxmin:idxmax) = 1;
[idxmin, idxmax] = get_box(400, 500, wavelengths);
Rideal(idxmin:idxmax) = 1;
[idxmin, idxmax] = get_box(500, 600, wavelengths);
Rideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(600, 700, wavelengths);
Rideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(700, 800, wavelengths);
Rideal(idxmin:idxmax) = 1;



Tideal = ones(size(wavelengths));
[idxmin, idxmax] = get_box(300, 400, wavelengths);
Tideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(400, 500, wavelengths);
Tideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(500, 600, wavelengths);
Tideal(idxmin:idxmax) = 1;
[idxmin, idxmax] = get_box(600, 700, wavelengths);
Tideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(700, 800, wavelengths);
Tideal(idxmin:idxmax) = 0;




% [idxmin, idxmax] = get_box(600, 700, wavelengths);
% Rideal(idxmin:idxmax) = 1;

R = @(w)get_R(w,ff,Escat, Escat0, Escat180);
A = @(w)get_A(w,ff,Eabs);
T = @(w)(1-get_R(w,ff,Escat, Escat0, Escat180)-get_A(w,ff,Eabs));

Ideal = Tideal;
Param = T;

f = @(w)norm(Ideal - Param(w)).^2;


w0 = ones(size(center_radiis))./length(center_radiis);
lb = zeros(size(w0));
ub = ones(size(w0));
Aeq = ones(size(w0));
beq = 1;

w = fmincon(f,w0,[],[],Aeq,beq,lb,ub);


Opt = Param(w);
norm(Ropt).^2;

%%

sizef = 20;


figure, 
subplot(2,3,1)
hold on 
bar(center_radiis, w0)
bar(center_radiis, w)
hold off
legend('Original', 'Optimized')
ylabel('Probability')
title('Area distribution')
lda = wavelengths(idxlda);
pbaspect([1 1 1])
box on 
set(gca, 'FontSize', sizef)


subplot(2,3,2)
plot(lda, 100.*T(w0), 'LineWidth', 6)
hold on 
plot(lda, 100.*T(w), 'LineWidth', 6)
plot(lda, 100.*Tideal, 'LineWidth', 6)
hold off
legend('Original','Optimized','Ideal')
xlim([min(lda), max(lda)])
ylabel('Transmission (%)') 
title('Film transmission')
pbaspect([1 1 1])
box on 
set(gca, 'FontSize', sizef)

subplot(2,3,3)
plot(lda, 100.*A(w0), 'LineWidth', 6)
hold on 
plot(lda, 100.*A(w), 'LineWidth', 6)
hold off
legend('Original','Optimized')
xlim([min(lda), max(lda)])
ylabel('Absorption (%)')
title('Film absorption')
pbaspect([1 1 1])
box on 
set(gca, 'FontSize', sizef)

subplot(2,3,4)
hold on 
for idx = 1:length(center_radiis)
    plot(lda, log10(FBR(:,idx)), 'LineWidth', 6)
end
hold off
xlim([min(lda), max(lda)])
ylabel('Log_{10} forward-to-backward ratio')
title('Isolated particle')
pbaspect([1 1 1])
box on 
set(gca, 'FontSize', sizef)

subplot(2,3,5)
hold on 
for idx = 1:length(center_radiis)
    plot(lda, Escat(:,idx), 'LineWidth', 6)
end
hold off
xlim([min(lda), max(lda)])
ylabel('Scattering efficiency')
title('Isolated particle')
pbaspect([1 1 1])
box on 
set(gca, 'FontSize', sizef)

subplot(2,3,6)
hold on 
for idx = 1:length(center_radiis)
    plot(lda, Eabs(:,idx), 'LineWidth', 6)
end
hold off
xlim([min(lda), max(lda)])
ylabel('Absorption efficiency')
title('Isolated particle')
pbaspect([1 1 1])
box on 
set(gca, 'FontSize', sizef)



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





































