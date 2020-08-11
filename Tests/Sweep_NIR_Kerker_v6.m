
clear
clc

core = '/home/parkerwray';
%addpath(genpath(strcat(core, '/hypnos/Codes/MSTMwrapper'))); % Keep these the same
addpath(genpath(strcat(core,'/hypnos/Codes/Matlab Functions')));  % Keep these the same
%addpath(genpath(strcat(core,'/hypnos/Codes/randomparticles')));  % Keep these the same
addpath(genpath(strcat(core, '/hypnos/Codes/Matlab_Mie_Core_Shell_Coefficients'))); % Keep these the same


x = 0.1:0.002:7;
m = [4];

nangsteps = 180;
opt = 1;
order = 10;
f = cell(length(x),1);


%% Generate the possible "reflection" profiles from an isolated particle

% Generate all the possible reflection profiles based on size parameter (x)
% and refractive index (m). The code assumes no loss, i.e., m = n+ik, k =0. 
% These reflection profiles form the basis at which you will try to
% construct an arbitrary reflection profile R = sum(Rj)

qa = zeros(length(x), length(m));
qs = qa;
FBR = qa;
a1 = qa;
a2 = qa;

tic;
for (idx2 = 1:1:length(m))
    parfor (idx = 1:1:length(x),5)
        f = pw_miecoated_Qs_FFs_Coeffs_v3(order, m(idx2), m(idx2), x(idx), x(idx), opt, nangsteps);    
        qa(idx,idx2) = f.qabs;
        qs(idx,idx2) = f.qsca;
        FBR(idx,idx2) = f.SL(1)/f.SL(90);
        If(idx, idx2) = f.SL(1);
        Ib(idx, idx2) = f.SL(90);
        a1(idx,idx2) = f.an(1);
        b1(idx,idx2) = f.bn(1);
        a2(idx,idx2) = f.an(2);
        b2(idx,idx2) = f.bn(2);
        a3(idx,idx2) = f.an(3);
        b3(idx,idx2) = f.bn(3);
        R(idx,idx2) = qs(idx, idx2)./(1+FBR(idx, idx2));
        A(idx, idx2) = qa(idx,idx2);
        T(idx,idx2) = qs(idx, idx2)-R(idx,idx2);
        
    end    
end
toc/60





%% Look at reflection function basis and identify region of interest 

% We care about basis functions such that: 
% R(lda) = sum(Rj(radii, lda, m)). 
% In reality we have Rj(x,m) where x = 2*pi*radii/lda.
% This allows or optimizer to pick the best radii (not just from a discrete
% set of pre-chosen radii) for a particular lda range. There are some
% subtlties in the implementation of this. First: for each radii, the data
% will span a different domain of wavelengths. Therefore, we need to
% simulate x larger than we actually desire. When converting to lda for a
% particular r, we are grabing a different subset of data within x. This
% smart truncation is performed in the optimization code. But, you need to
% feed the code a range of x where it can properly perform the subset
% truncation. Therefore, you should visibly check the data. Make sure that
% the reflections 

% Find radii that for upper and lower limit based on wavelength range
clc
lda = [2,6].*1000;
%r = outer_product(x,lda)./(2*pi);

ldamin = min(lda);
ldamax = max(lda);
[~, idx_dipole] = min(abs(0.5-x));
[~, idx_multipole] = min(abs(3-x));

plot_flag = 1;


for idxm = 1:length(m)
    
    % Find location of reflection maximum within the region between 
    % the Rayleigh and Multipole limit.
    [vmax, idx_x_max] = max(R(idx_dipole:idx_multipole,idxm));
    
    % Find location of reflection minimum within the region between 
    % the Rayleigh and Multipole limit.    
    [vmin, idx_x_min] = min(R(idx_dipole:idx_multipole,idxm));
    
    % Get the size parameters associated with reflection minimum and
    % maximum
    xmax(idxm) = x(idx_x_max+idx_dipole);
    xmin(idxm) = x(idx_x_min+idx_dipole);
    
    % Get the radius range such that you can choose from a range of
    % particles that have reflection minimums and maximums spanning the
    % space
    rmin(idxm) = xmin(idxm).*ldamin./(2*pi); % We want to capture the peak at the shortes wavelength
    rmax(idxm) = xmax(idxm).*ldamax./(2*pi); % we want to capture the backward at the longest wavelength

end

plot_flag = 1;
if plot_flag ==1
    figure,
    for idxm = 1:length(m)   
        subplot(length(m),1,idxm)
        hold on 
        plot(x, R(:,idxm));
        plot(x(idx_dipole:idx_multipole), R(idx_dipole:idx_multipole,idxm));
        xline(xmin(idxm));
        yline(vmin);
        xline(xmax(idxm));
        yline(vmax);
        hold off
        legend('Full sim. region','Min/max search region')
        ylabel('Pseudo reflection (a.u.)')
        xlabel('Size parameter (x)')  
    end
    figure,
    for idxm = 1:length(m)   
        subplot(length(m),1,idxm)
        hold on 
        for r = linspace(rmin(idxm), rmax(idxm), 5)  %linspace(150,295,2)
            lda = 2*pi*r./x;
            plot(lda/1000,R(:,idxm));
        end
        hold off
        ylabel('Pseudo reflection (a.u.)')
        xlabel('Wavelength (um)') 
        xlim([ldamin, ldamax]/1000)
    end   
end

 

%% Make the refleciton profile that you want to generate


idxm = 1;
ff = 0.3;

figure, hold on 
r = linspace(rmin(idxm), rmax(idxm), 50);
lda = 2*pi*r(15)./x;
plot(lda/1000, ff.*0.5.*R(:,idxm));
lda = 2*pi*r(45)./x;
plot(lda/1000, ff.*0.5.*R(:,idxm));

lda = linspace(ldamin, ldamax, 300);
% Make desired reflection distribution
Rideal = zeros(size(lda)).';
[idxmin, idxmax] = get_box(2000, 3000, lda);
Rideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(3000, 4000, lda);
Rideal(idxmin:idxmax) = 0;
[idxmin, idxmax] = get_box(4000, 5000, lda);
Rideal(idxmin:idxmax) = 0.8;
[idxmin, idxmax] = get_box(5000, 6000, lda);
Rideal(idxmin:idxmax) = 0.8; 
Rideal = Rideal.*1;

plot(lda/1000, Rideal);


r1 = r(15);
r2 = r(45);


%
clc
% Let v be a vector of elements [ff, radiis, radii weights]
% Where radiis = nx1, and radii weights = nx1
% For fill fraction -> 0.2 < ff < 0.4
% For all radiis -> r lower bound < r < r upper bound
% Note: the range of radii will depend on refractive index per code above.
% For all weights -> 0 < w < 1 and sum(w) = 1
clearvars ffmin ffmax resolution
clearvars minimum_radii maximum_radii
clearvars minimum_weights maximum_weights
clearvars lb ub Aeq beq v0
% Define constraints
ffmin = 0.1;
ffmax = 0.4;
Lparticles = 5; %5;
Lm = length(m);

minimum_radii = [];
maximum_radii = [];
for idxm = 1:length(m) 
    minimum_radii = [minimum_radii; rmin(idxm).*ones(Lparticles,1)];
    maximum_radii = [maximum_radii; rmax(idxm).*ones(Lparticles,1)];
end

minimum_weights = zeros(size(minimum_radii));
maximum_weights = ones(size(maximum_radii));

% Define upper and lower bound for x
%lb = [ffmin; minimum_radii; minimum_weights];
%ub = [ffmax; maximum_radii; maximum_weights];
lb = [ffmin; minimum_weights; minimum_weights];
ub = [ffmax; maximum_weights; maximum_weights];

% We want Aeq*x = sum(w) = beq = 1
Aeq = [0; zeros(size(minimum_radii)); ones(size(minimum_weights))].';
beq = 1;

ff0 = ffmin+(ffmax-ffmin).*rand(1);

r0 = minimum_radii + (maximum_radii-minimum_radii).*rand(Lparticles*Lm,1); 
r0 = (r0-minimum_radii)./(maximum_radii-minimum_radii);

w0 = rand(Lparticles*Lm,1);
w0 = w0./sum(w0);

v0 = [ff0; r0; w0];

% r0 = minimum_radii;
% r0(7) = r1;
% r0(8) = r2;
% r0 = (r0-minimum_radii)./(maximum_radii-minimum_radii);
% v0 = [0.3; r0; 0;0;0; 0;0;0; 0.3;0.7;0];


% v0 = [(ffmin+ffmax)./2; ...
%       (minimum_radii+maximum_radii)./2;...
%       ones(size(minimum_weights))./length(minimum_weights)];
% 
% v0 = [(ffmin+ffmax)./2; ...
%       (minimum_radii);...
%       ones(size(minimum_weights))./length(minimum_weights)];  
%   
  
size_m = length(m);

lda = linspace(ldamin, ldamax, 300);

[qs_reshaped, If_reshaped, Ib_reshaped] = ...
    reshape_for_Lparticles(Lparticles, Lm, qs, If, Ib);


g = @(v)get_R_optimizer(v, ...
                    Lparticles,...
                    Lm,...
                    lda,...
                    x,...
                    qs_reshaped,...
                    If_reshaped,...
                    Ib_reshaped,...
                    maximum_radii,...
                    minimum_radii);
R0 = g(v0);
% figure, hold on 
plot(lda/1000, R0,':')
% plot(lda, Rideal)
hold off
xlim([2,6])

%
clc
%
%v0 = vopt;
% Optimizer options
options = optimoptions('fmincon');
options.Algorithm = 'interior-point';
options.StepTolerance = 1e-20;
options.ConstraintTolerance = 1e-20;
options.OptimalityTolerance = 1e-20;
options.Display = 'iter-detailed';
options.MaxIterations = 100000;
options.MaxFunctionEvaluations = 10000;
% options.FiniteDifferenceStepSize = (ub-lb)/10;
% options.DiffMinChange = 0;

tic

                
f = @(v)(norm(Rideal-g(v)).^2);
%f = @(v)sum(abs(Rideal-g(v)));                                
[vopt, vval] = fmincon(f,v0,[],[],Aeq,beq,lb,ub,[], options);
toc/60

% Get results of the optimal parameters and plot them vs. the seed. 
%clc
%plot_v_opt_vs_seed(Lparticles, Lm, vopt, v0);


Ropt = g(vopt);

[ffopt, ropt, wopt] = get_parameters(vopt, Lparticles, Lm, minimum_radii, maximum_radii);

%%
saveflag = 0;
core = 'Custom2';

fig = figure; 
[rs, Is] = sort(ropt(:,1));
bar(wopt(Is,1).*100)
xticklabels(num2str(round(rs)))
ylabel('Probability (%)')
xlabel('Radius (nm)')
title(['m = ', num2str(m(1)),'+j0'])
set(gca,'FontSize',30)
pbaspect([1,1,1])
box on 
set(gcf, 'Position', get(0, 'Screensize'));
if saveflag
savefig([core,'_Dist1.fig'])
saveas(gcf,[core,'_Dist1.tif']) 
end

if length(m)>1
figure,
[rs, Is] = sort(ropt(:,2));
bar(wopt(Is,2).*100)
xticklabels(num2str(round(rs)))
ylabel('Probability (%)')
xlabel('Radius (nm)')
title(['m = ', num2str(m(2)),'+j0'])
set(gca,'FontSize',30)
pbaspect([1,1,1])
box on 
set(gcf, 'Position', get(0, 'Screensize'));
if saveflag
savefig([core,'_Dist2.fig'])
saveas(gcf,[core,'_Dist2.tif']) 
end
if length(m)>2
figure,
[rs, Is] = sort(ropt(:,3));
bar(wopt(Is,3).*100)
xticklabels(num2str(round(rs)))
ylabel('Probability (%)')
xlabel('Radius (nm)')
title(['m = ', num2str(m(3)),'+j0'])
set(gca,'FontSize',30)
pbaspect([1,1,1])
box on 
set(gcf, 'Position', get(0, 'Screensize'));
if saveflag
savefig([core,'_Dist3.fig'])
saveas(gcf,[core,'_Dist3.tif']) 
end
end
end


figure, 
hold on 
plot(lda/1000, 100.*Ropt,'k-','LineWidth',8)
plot(lda/1000, 100.*R0,'b-.','LineWidth',8)
plot(lda/1000, 100.*Rideal,'r:','LineWidth',8)
hold off
xlabel('Wavelength (um)')
ylabel('Reflection (%)')
legend(['Optimized',newline, '(ff = ',num2str(round(ffopt*100)),'%)'],...
    ['Seed',newline, '(ff = ', num2str(round(ff0*100)), '%)'],...
    'Target profile')
legend boxoff
set(gca,'FontSize',30)
pbaspect([1,1,1])
box on 
set(gcf, 'Position', get(0, 'Screensize'));
if saveflag
savefig([core,'_Reflection.fig'])
saveas(gcf,[core,'_Reflection.tif']) 
end
%%
function r = scale_r(r0, min_r, max_r)

r = r0.*(max_r-min_r)+min_r;

end

function [ff, r, w] = get_parameters(v, Lparticles, Lm, minimum_radii, maximum_radii)

ff = v(1);
rl = v(2:Lparticles*Lm+1);
wl = v(Lparticles*Lm+2:end);

rl = scale_r(rl, minimum_radii, maximum_radii);

for idx = 1:length(rl)
    [idx_particle, idx_m] = ind2sub([Lparticles, Lm], idx);
    r(idx_particle, idx_m) = rl(idx);
    w(idx_particle, idx_m) = wl(idx);
end



end


function [qs, If, Ib] = reshape_for_Lparticles(Lparticles, Lm, qs, If, Ib)

    % Repeat the qs and FBR data based on the number of distinct radii you
    % want to consider (i.e., the "resolution" of the probability
    % distribution). In this case, the qs and fbr data are changing with
    % refractive index (m = n+ik). For each distinct refractive index we 
    % will consider Nparticles of variable radii, which the optimizer will
    % find for us. This loop says, for each m, repeat the qs and fbr data
    % for Nparticles.

%     rfactor = ((maximum_radii-minimum_radii)+minimum_radii).';
    
    qsdummy = [];
    Ifdummy = [];
    Ibdummy = [];
%     rfactordummy = [];
    for idx = 1:Lm
        qsdummy = [qsdummy,repmat(qs(:,idx),1,Lparticles)];
        Ifdummy = [Ifdummy,repmat(If(:,idx),1,Lparticles)];
        Ibdummy = [Ibdummy,repmat(Ib(:,idx),1,Lparticles)];
%         rfactordummy = [rfactordummy,repmat(rfactor(:,idx),1,Lparticles)];
    end
    clearvars qs If Ib
    qs = qsdummy.';
    If = Ifdummy.';
    Ib = Ibdummy.';

    % qs and fbr are now size [m*Nparticles, x]
   
end


function [qs, If, Ib] = get_domain_for_different_radiis(radiis, lda, x, qs, If, Ib)
    
    % minimum and maximum lda range you want to consider in your simulation
    % as well as the lda spacing you want to consider. (i.e., the data all
    % has the same x spacing. But that is not the same as having the same
    % lda spacing. You want to interpolate to make the same lda spacing). 
    ldamin = min(lda);
    ldamax = max(lda);
    Llda = length(lda);
    
    % lda is dim radiis x size param
    ldas = 2*pi*radiis./x;
    % For each new radii there will be a different range of lda that is
    % spanned because x is static. We dont need all of this data. We only
    % need the data that corresponds to the wavelength range we are
    % intersted in. 
    
    [~, idx_ldamin] = min(abs(ldas-ldamin),[],2);
    [~, idx_ldamax] = min(abs(ldas-ldamax),[],2);
    
    idx_xmin = idx_ldamax; %[m*Nparticles x 1]
    idx_xmax = idx_ldamin; %[m*Nparticles x 1]

    % These idices determine the index values in x we want to grab data
    % from for this particular radii within our wavelength range of
    % interest. I.e., for each radii this tells us the x data range we want
    % to grab.
   
    
    % For each radii, grab the correct domain of x data and interpolate it
    % such that all data is the same size. The resulting matrix will be 
    % [Number desired lda pts within desired window x number of index pts.]
    for idx = 1:size(qs,1)
       
       % Grab the relevant lda domain for particle idx.
       domain = idx_xmin(idx):idx_xmax(idx);
       
       qs_dummy(idx,:) = interp1(ldas(idx, domain),...
                        qs(idx,domain),...
                        lda, 'spline','extrap'); 
       
       If_dummy(idx,:) = interp1(ldas(idx, domain),...
                        If(idx,domain),...
                        lda, 'spline','extrap'); 
                    
       Ib_dummy(idx,:) = interp1(ldas(idx, domain),...
                        Ib(idx,domain),...
                        lda, 'spline','extrap'); 
      % Plot values for debugging/checking purposes.              
%       figure, 
%       plot(ldas(idx, domain),qs(idx,domain))
%       hold on 
%       plot(lda, qs_dummy,':', 'LineWidth', 6)
%       hold off              
%       figure, 
%       plot(ldas(idx, domain),fbr(idx,domain))
%       hold on 
%       plot(lda, fbr_dummy,':', 'LineWidth', 6)
%       hold off              
            
    end                       
    clearvars qs If Ib
    qs = qs_dummy;
    If = If_dummy;
    Ib = Ib_dummy;
end


function R = get_R_optimizer(v, Lparticles, Lm, lda, x, qs, If, Ib, maximum_radii,minimum_radii)
    
    % x = [ff, radiis, radii weights]
    % R = ff*qs/(1+FBR)
    % qs and FBR are given as [size param, m]

    ff = v(1);
    radiis = (v(2:Lparticles*Lm+1)).*(maximum_radii-minimum_radii)+minimum_radii;
    weights = (v((Lparticles*Lm)+2:end)).';
       
    [qs, If, Ib] = get_domain_for_different_radiis(radiis, lda, x, qs, If, Ib); 


    wqs = (weights*qs).';
    wIf = (weights*If).';
    wIb = (weights*Ib).';
    wfbr = wIf./wIb;
    
    R = (ff.*(wqs)./(1+wfbr));
   
    R(R>1) = 1;
    R(R<0) = 0;
end






function [idxmin, idxmax] = get_box(ldamin, ldamax, wavelengths)


[~,idxmin] = min(abs(wavelengths-ldamin));
[~,idxmax] = min(abs(wavelengths-ldamax));

end



function plot_v_opt_vs_seed(Lparticles, Lm, vopt, v0)
ffopt = vopt(1);
ropt = vopt(2:Lparticles*Lm+1);
wopt = vopt(Lparticles*Lm+2:end);
ff0 = v0(1);
r0 = v0(2:Lparticles*Lm+1);
w0 = v0(Lparticles*Lm+2:end);
100*(vopt-v0)./v0;

figure, 
subplot(2,1,1)
yyaxis left
hold on 
plot(ropt)
plot(r0)
hold off
ylabel('Particle radius') 
yyaxis right
plot(100*(ropt-r0)./r0)
ylabel('Fractional change from seed')
xlabel('Particle number')
legend('Optimized','Seed','% change')
subplot(2,1,2)
yyaxis left
hold on 
plot(wopt)
plot(w0)
hold off
ylabel('Particle weight') 
yyaxis right
plot(100*(wopt-w0)./w0)
ylabel('Fractional change from seed')
xlabel('Particle number')
legend('Optimized','Seed','% change')

disp(['Optimized fill fraction is ', num2str(100*ffopt),'%'])
disp(['Seed fill fraction was ', num2str(100*ff0),'%'])
disp(['Fractional change of ', num2str(100*(ffopt-ff0)./ff0),'%'])    

end




