
% clearvars -except File_out File_scat
% clc

%% Get Files and data 
if ~exist('File_out')
    [File_out, FileName, PathName] = get_file();
end
if ~exist('File_scat')
[   File_scat, FileName, PathName] = get_file(PathName);
end


[cluster, sphere, excitation] = export_output_file(File_out);
[sphere_coeffs, alpha, beta, cbeam, m_med, Nspheres, Nequs] = export_mstm_scattering_coeffs(File_scat);


%%
theta = [cluster.S.theta];
I_unpolarized = [cluster.S.s11];
I_parallel = [cluster.S.s11]+[cluster.S.s12];
I_perpendicular = [cluster.S.s11]-[cluster.S.s12];

figure, 
plot(theta, log10(I_unpolarized), 'LineWidth', 6)
hold on 
plot(theta, log10(I_parallel), 'LineWidth', 6)
plot(theta, log10(I_perpendicular), 'LineWidth', 6)
hold off
xlabel('Theta (deg)')
ylabel('Log_{10} Far field intensity')
legend('Unpolarized','Parallel','Perpendicular')
xlim([0,180])
pbaspect([1 1 1])
set(gca, 'FontSize',24)

%%
figure, 
polar(theta*pi/180, I_parallel)
hold on 
polar(theta*pi/180, I_perpendicular)
hold off


%%


%%% ******* Use full!!!
% 
% a1_s = sphere_coeffs.a_te{1};
% a2_s = sphere_coeffs.a_tm{1};
% 
% 
% 
% a1 = 0.5.*(a1_s+a2_s);
% a2 = 0.5.*(a1_s-a2_s);


% a1 = 2.*a1;
% b1 = 2.*a2;



%%

% b_tilde = sphere_coeffs.a_te;
% a_tilde = sphere_coeffs.a_tm;
% 
% k = 0.0137;
% 
% 
% 
% 
% 
% %
% a1_tilde = a_tilde{1};
% b1_tilde = b_tilde{1};
% % 
% % 
% a_test_1 = sum(a1_tilde);
% a_test_2 = sum(b1_tilde);
% 
% a1 = sum(0.5.*(a1_tilde+b1_tilde));
% b1 = sum(0.5.*(a1_tilde-b1_tilde));
% 
% 
% g_tilde = sphere_coeffs.f_te;
% f_tilde = sphere_coeffs.f_tm;
% 
% g1_tilde = g_tilde{1};
% f1_tilde = f_tilde{1};
% 
% g1 = sum(0.5.*(g1_tilde+f1_tilde));
% f1 = sum(0.5.*(g1_tilde-f1_tilde));




