function spheres = Dual_Sphere_Anisotropic_Kerker_20_06_23(wavelengths)

r_lens = 10*wavelengths/2*pi;
r_abs = 0.15*r_lens;

spheres = [r_lens, 0, 0, 0, 2, 0;...
            r_abs, 0, 0, r_lens-r_abs-0.001, 3, 4];

end
