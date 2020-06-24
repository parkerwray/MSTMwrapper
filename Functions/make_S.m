function [Ipar, Iper, S1, S2, S3, S4] = make_S(A, B, Ax, Bx, theta)

% Structures A, B, Ax, Bx, have the dimensions:
% [wavelength, polarization, order, mode]
% If your structure has the dimensions:
% [simulations, wavelength, polarization, order, mode]
% you need to pass this function for each simulation.
if length(size(A)) > 4
    A = squeeze(A);
    B = squeeze(B);
    Ax = squeeze(Ax);
    Bx = squeeze(Bx);
end
size_theta = length(theta);

[size_lda, size_pol, size_n, size_m] = size(A);
S1 = zeros(size_lda,size_theta);
S2 = S1;
S3 = S1;
S4 = S1;

par = 1;
per = 2;
for n = 1:size_n
    [tau_n, pi_n] = make_legendre(n,theta);

    for m = 1:n+1
        S1 = S1+(A(:,per,n,m).*pi_n(m,:)+B(:,per,n,m).*tau_n(m,:));
        S2 = S2+(A(:,par,n,m).*tau_n(m,:)+B(:,par,n,m).*pi_n(m,:));
        S3 = S3+(Ax(:,per,n,m).*tau_n(m,:)+Bx(:,per,n,m).*pi_n(m,:));
        S4 = S4+(Ax(:,par,n,m).*pi_n(m,:)+Bx(:,par,n,m).*tau_n(m,:));
    end
end

[Ipar, Iper] = make_I(S1, S2, S3, S4);


function [Ipar, Iper] = make_I(S1, S2, S3, S4)

    Iper = NaN(size(S1));
    Ipar = Iper;
    for idx = 1:size(S1,1)
        for idx2 = 1:size(S1,2)
            Iper(idx,idx2) =  abs(S2(idx,idx2)).^2+...
                               abs(S3(idx,idx2)).^2+...
                               dot(S2(idx,idx2),conj(S3(idx,idx2)))+...
                               dot(conj(S2(idx,idx2)),(S3(idx,idx2)));

            Ipar(idx,idx2) =  abs(S4(idx,idx2)).^2+...
                              abs(S1(idx,idx2)).^2+...
                              dot(S4(idx,idx2),conj(S1(idx,idx2)))+...
                              dot(conj(S4(idx,idx2)),(S1(idx,idx2)));
        end
    end
    


end



end


