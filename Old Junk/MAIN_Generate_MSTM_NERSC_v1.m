
parentdir = uigetdir('$ME','Specify Save Directory');


%% SIMULATION INPUT PARAMETERS

beam_type = 0; % 0 = plane wave, 1 = Gaussian beam

if beam_type == 1
    beam_waist = ceil(1./(0.19.*k)); % Spot RADIUS at focus
    Zr = (pi.*beam_waist.^2)./lda; % Rayleigh range in [nm]
    disp(['The input beam diameter is ', num2str(2.*beam_waist), 'nm at focus'])
    disp(['The depth of field is ', num2str(2.*Zr), 'nm'])
    disp(['The unitless parameter is ', num2str(1./(k.*beam_waist)), ' (should be <= 0.2)'])
else
    beam_waist = 0;
end

%beam_width = lda;   % 0 = plane wave, 1/w*k = Gaussian beam, beam_width = width [nm];
alpha = 0; % input wave polar angle
beta = 0; % input wave venith angle
pol = 0; % polarization state
near_field_cords = [-3000,-3000,3000,6000]; % [-x ,-y ,x ,y] Coordinates [unit nm]
near_field_resolution = 0.5;
wavelengths = 300:10:400;

medium_m = 1;
spheres = [500, 0, 0, 0, 5, 3]; 
Nspheres = 1;
for idx = 1:length(wavelengths)
fname= strcat('mstm',... %GIVE A FILE NAME!
        '_lda_', num2str(wavelengths(idx)),'nm');
k = 2*pi/wavelengths(idx);
make_mstm_NERSC_job(parentdir, fname, sphere, Nspheres, medium_m, k,...
    beam_type, beam_waist,alpha, beta, pol, near_field_cords,...
    near_field_resolution);
end
  