%% Move new data to secondary variable
A2 = A;
Ax2 = Ax;
B2 = B;
Bx2 = Bx;
Ipar2 = Ipar;
Iper2 = Iper;
qa2 = qa;
qe2 = qe;
qsd2 = qsd;
qsi2 = qsi;

%% Clear variables just to make sure import works
clearvars A Ax B Bx Ipar Iper qa qe qsd qsi

%% Concatenate raw data
dim_N_modes = 2;
dim_N_other = 3;

A  = cat(dim_N_modes, A , A2 );
Ax = cat(dim_N_modes, Ax, Ax2);
B  = cat(dim_N_modes, B , B2 );
Bx = cat(dim_N_modes, Bx, Bx2);

Iper = cat(dim_N_other, Iper, Iper2);
Ipar = cat(dim_N_other, Ipar, Ipar2);

qa  = cat(dim_N_other, qa , qa2 );
qe  = cat(dim_N_other, qe , qe2 );
qsd = cat(dim_N_other, qsd, qsd2);
qsi = cat(dim_N_other, qsi, qsi2);

%% Recreate statistics
dim_wavelengths = 3;

[Astats] = get_mode_statistics(A, dim_N_modes);
[Bstats] = get_mode_statistics(B, dim_N_modes);
[Axstats] = get_mode_statistics(Ax, dim_N_modes);
[Bxstats] = get_mode_statistics(Bx, dim_N_modes);

Astats.energy = get_range(Astats.mean.mag,dim_wavelengths);
Bstats.energy = get_range(Bstats.mean.mag,dim_wavelengths);
Axstats.energy = get_range(Axstats.mean.mag,dim_wavelengths);
Bxstats.energy = get_range(Bxstats.mean.mag,dim_wavelengths);


[qestats] = get_statistics(qe, dim_N_other);
[qastats] = get_statistics(qa, dim_N_other);
[qsistats] = get_statistics(qsi, dim_N_other);
[qsdstats] = get_statistics(qsd, dim_N_other);

%%


I = mean(Ipar(:,:,:,:,:),3)+mean(Iper(:,:,:,:,:),3)./2;

[qestats1]  = get_statistics(qe (:,:,1:25,:), dim_N_other);
[qastats1]  = get_statistics(qa (:,:,1:25,:), dim_N_other);
[qsistats1] = get_statistics(qsi(:,:,1:25,:), dim_N_other);
[qsdstats1] = get_statistics(qsd(:,:,1:25,:), dim_N_other);
I1 = mean(Ipar(:,:,1:25,:,:),3)+mean(Iper(:,:,1:25,:,:),3)./2;


[qestats2]  = get_statistics(qe (:,:,26:50,:), dim_N_other);
[qastats2]  = get_statistics(qa (:,:,26:50,:), dim_N_other);
[qsistats2] = get_statistics(qsi(:,:,26:50,:), dim_N_other);
[qsdstats2] = get_statistics(qsd(:,:,26:50,:), dim_N_other);
I2 = mean(Ipar(:,:,26:50,:,:),3)+mean(Iper(:,:,26:50,:,:),3)./2;

%% SAVE!!!
core = '/home/elipaul';
savedir = uigetdir(core,' Specify Save Directory'); % for Linux
old_location = cd(savedir);
mkdir('saved_data');
cd('saved_data');
save('modes.mat','A','B','Ax','Bx', '-v7.3');
save('Is.mat', 'Ipar', 'Iper', '-v7.3');
save('qs.mat', 'qe','qa','qsi','qsd','-v7.3');
%save('sphere_distributions.mat', '-v7.3');
save('simulation_parameters.mat', 'fill_fractions', 'center_radiis', 'Ndistributions', 'max_order', 'wavelengths', '-v7.3');    % UPDATE THIS LINE WITH THE RELEVANT PARAMETERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
save('modes_stats.mat','Astats','Bstats','Axstats','Bxstats','-v7.3');
save('qs_stats.mat','qestats','qastats','qsistats','qsdstats','-v7.3');
save('wopt.mat', 'wopt', '-v7.3');

