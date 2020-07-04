

function [A, B, Ax, Bx, size_param,...
            Apar, Bpar, Axpar, Bxpar,...
            Aper, Bper, Axper, Bxper] = ...
                process_convert_modes(sphere_modes_par, sphere_modes_per)
%{ 
This function converts MSTM modes to Apnm, Bpnm, Axpnm, Bxpnm, where 
p is size 2, [1 = parallel, 2 = perpendicular]
n is size Ncritical, [1, 2, ...,  Ncritical]
m is size Ncritical+1, [1 -> m = 0, 2 -> m = 1, ...]

The code scales the modes such that they maintain the conventions shown
in Mie theory. I.e., the modes are scaled by sqrt(2/(2n+1)). This creates a
perfect match when comparing modes from Mie theory and MSTM code. 

Test of comparing the modes and efficiencies for a single sphere using MSTM
and Mie theory code show a near-perfect match. Dimer tests show that the
scaled MSTM modes for the particle at the origin will continue to reproduce 
the scattering cross section of the particle at the origin. Therefore, the
scalling works for coupled particles as well. 
%}


[len_lda, len_sims] = size(sphere_modes_par);
if len_lda == 1 && len_sims > 1
    % This funciton assumes [n 1] or [n m]
    % if sphere_modes is [1 n] flip it.
    sphere_modes_par = sphere_modes_par.';
    sphere_modes_per = sphere_modes_per.';
    dummy = len_lda;
    len_lda = len_sims;
    len_sims = dummy;
    clearvars dummy
end

Max_Order = 10;    

Apar= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Bpar= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Axpar= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Bxpar= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Aper= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Bper= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Axper= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 
Bxper= zeros(len_sims, len_lda, Max_Order, Max_Order+1); 

% Get Parallel Modes
    for idx_sims = 1:len_sims
        for idx_lda = 1:len_lda
            Nmax = sphere_modes_par{idx_lda,idx_sims}.Nsphere_order;
            for idx_n = 1:Nmax
            polarization = 0; % Flag for parallel is 0
            [A,B,Ax,Bx] = ...
                unpack_modes(sphere_modes_par{idx_lda,idx_sims},idx_n,...
                                polarization);    
                for idx_m = 1:length(A)     
                    if isempty(A(idx_m))
                        Apar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Apar(idx_sims, idx_lda, idx_n, idx_m) = A(idx_m);
                    end  
                    if isempty(B(idx_m))
                        Bpar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bpar(idx_sims, idx_lda, idx_n, idx_m) = B(idx_m);
                    end                      
                    if isempty(Ax(idx_m))
                        Axpar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Axpar(idx_sims, idx_lda, idx_n, idx_m) = Ax(idx_m);
                    end  
                    if isempty(Bx(idx_m))
                        Bxpar(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bxpar(idx_sims, idx_lda, idx_n, idx_m) = Bx(idx_m);
                    end     
                    
                   
                end 
            end
        end
    end
% Get Perpendicular Modes
    for idx_sims = 1:len_sims
        for idx_lda = 1:len_lda
            Nmax = sphere_modes_par{idx_lda,idx_sims}.Nsphere_order;
            
            % Old code used "ka" for size parameter instead of "x". 
            % Account for both styles with a logic statement. 

            if isfield(sphere_modes_par{idx_lda,idx_sims}(1),'x')
                if ~isempty(sphere_modes_par{idx_lda,idx_sims}(1).x)
                size_param(idx_lda) =...
                    sphere_modes_par{idx_lda,idx_sims}(1).x;
                end
            end
            
            if isfield(sphere_modes_par{idx_lda,idx_sims}(1),'ka')
                if ~isempty(sphere_modes_par{idx_lda,idx_sims}(1).ka)
                size_param(idx_lda) =...
                    sphere_modes_par{idx_lda,idx_sims}(1).ka;
                end
            end

            
            for idx_n = 1:Nmax
            polarization = 1; % Flag for perpendicular is 1
            [A,B,Ax,Bx] = ...
                unpack_modes(sphere_modes_per{idx_lda,idx_sims},idx_n,...
                                polarization);
                for idx_m = 1:length(A)     
                    if isempty(A(idx_m))
                        Aper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Aper(idx_sims, idx_lda, idx_n, idx_m) = A(idx_m);
                    end  
                    if isempty(B(idx_m))
                        Bper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bper(idx_sims, idx_lda, idx_n, idx_m) = B(idx_m);
                    end                      
                    if isempty(Ax(idx_m))
                        Axper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Axper(idx_sims, idx_lda, idx_n, idx_m) = Ax(idx_m);
                    end  
                    if isempty(Bx(idx_m))
                        Bxper(idx_sims, idx_lda, idx_n, idx_m) = NaN;
                    else
                        Bxper(idx_sims, idx_lda, idx_n, idx_m) = Bx(idx_m);
                    end                        
                    
                end
            end
        end
    end
    
    A = cat(1, Apar, Aper);
    B = cat(1, Bpar, Bper);
    Ax = cat(1, Axpar, Axper);
    Bx = cat(1, Bxpar, Bxper);
    A = permute(A,[2,1,3,4]);
    B = permute(B,[2,1,3,4]);
    Ax = permute(Ax,[2,1,3,4]);
    Bx = permute(Bx,[2,1,3,4]);
 
 
 
function [A, B, Ax, Bx] = unpack_modes(modes, n, polarization)
% Parallel polarization = 0
% Perpendicular polarization = 1
A = [];
B = [];
Ax = [];
Bx = [];

        orig = modes(1); % Grab the sphere at the origin data
        modes_a_par = orig.a_te{n};
        modes_b_par = orig.a_tm{n};

        %Sort Coeffs
        m = (-n:1:n);
        A_neg_m_n = [];
        B_neg_m_n = [];
        A_m_n = [];
        B_m_n = [];
        for idx_dummy = 1:length(m)
            if m(idx_dummy)<0
                A_neg_m_n = [A_neg_m_n,  modes_a_par(idx_dummy)];
                B_neg_m_n = [B_neg_m_n,  modes_b_par(idx_dummy)];
            elseif m(idx_dummy)>0
                A_m_n = [A_m_n, modes_a_par(idx_dummy)];
                B_m_n = [B_m_n, modes_b_par(idx_dummy)];
            elseif m(idx_dummy) == 0
                A_0_n = modes_a_par(idx_dummy);
                B_0_n = modes_b_par(idx_dummy);
            end
        end
        
        % m = -n to n.  Therefore, the array for the negative values goes
        % [-4, -3, -2, -1] and the positive values goes [1, 2, 3, 4]. To 
        % add the right elements, we need the vectors to be in the right 
        % order. 
        A_neg_m_n = fliplr(A_neg_m_n);
        B_neg_m_n = fliplr(B_neg_m_n);
        
        % Scaling factor to make MSTM modes similar to the ones in Mie
        % theory. 
        scale_factor = sqrt(2/(2*n+1));
        
        % Initialize the first element with An0 and Bn0 for the 
        % cross-polarization, there is no m = 0, so this element is NAN. 
        A = A_0_n;
        B = B_0_n;
        Ax = 0; %NaN;
        Bx = 0; %NaN;
        
        for m = 1:n
          
            dummy =(((-1).^(n)).*(1i).^(n+1+polarization)).*...
                (A_m_n(m)+((-1).^(m+polarization)).*A_neg_m_n(m))...
                .*scale_factor; %.*

            A =...
                [A, dummy];

            dummy = ...              
                (((-1).^(n)).*(1i).^(n+1+polarization))...
                .*(B_m_n(m)+((-1).^(m+1+polarization)).*B_neg_m_n(m))...
                .*scale_factor;

            B =...
                [B, dummy];

            dummy = ...
                (((-1).^(n)).*(1i).^(n+1+polarization))...
                .*(A_m_n(m)+((-1).^(m+1+polarization)).*A_neg_m_n(m))...
                .*scale_factor;

            Ax =...
                [Ax, dummy];

            dummy = ...
                (((-1).^(n)).*(1i).^(n+1+polarization))...
                .*(B_m_n(m)+((-1).^(m+polarization)).*B_neg_m_n(m))...
                .*scale_factor;

            Bx =...
                [Bx, dummy];                   
        end
end







end










