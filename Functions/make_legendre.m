function[tau_nm, pi_nm]= make_legendre(n,theta)



% n = 10;
m = (0:1:n).';
% theta = (0:0.01:360)*pi/180;
idx_pi = find(theta == pi);
theta(theta == 0) = 0.0000001;
theta(theta == pi) = pi-0.0000001;
theta(theta == 2*pi) = 2*pi-0.0000001;
x = cos(theta);

P = legendre(n,x);
[size_m, size_x] = size(P);
Pprev = [legendre(n-1,x); zeros(1,size_x)];
nX = n.*repmat(x,size_m,1);
M = repmat(m, 1, size_x);
npM = n.*ones(size(M))+M;


tau_nm = sin(theta).*(nX.*P-npM.*Pprev)./(x.^2-1); 
pi_nm = M.*P./sin(theta);

tau_nm(isnan(tau_nm)) = 0;
pi_nm(isnan(pi_nm)) = 0;



for idx_m = 1:length(m)
%     figure, 
%     hold on 
%     plot(tau_nm(idx_m,:))
    if ~mod(idx_m,2)
        tau_nm(idx_m,theta>pi) = -tau_nm(idx_m,theta>pi);
        pi_nm(idx_m,theta<pi) = -pi_nm(idx_m,theta<pi);
    end
%     plot(tau_nm(idx_m,:),':','LineWidth',4)
%     hold off
%     title(['m = ', num2str(idx_m-1)])
%     legend('pre-correction', 'post-correction')
    scale_n = ((2*n+1)/(n*(n+1))); %sqrt((n*(n+1))./(2*n+1)); sqrt((2*n+1)/2); sqrt((n*(n+1))./(2*n+1));  %2.125; %.*sqrt((2)/(2*n+1));
    scale_m = 1; %sqrt(factorial(n-(m-1))/factorial(n+(m-1)));
    scale_nm= scale_n.*scale_m;

    tau_nm(idx_m,:)=(tau_nm(idx_m,:)).*scale_nm;
    pi_nm(idx_m,:)=(pi_nm(idx_m,:)).*scale_nm;
end
end
