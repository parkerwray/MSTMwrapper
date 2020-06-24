


clc;
k = 21;
r = 44;


sphere_cords = (sphere_coeffs_par(k,:));


%%


for j = 1:length(sphere_cords)
    cords = [];
    cords(1,:) = [sphere_cords{j}.loc_x];
    cords(2,:) = [sphere_cords{j}.loc_y];
    cords = cords.';
    lcords = size(cords,1);
    
    for i = 1:lcords

        P = cords;
        P(i,:) = [];
        PQ = cords(i,:);

         [idx(i), dist(i)] = dsearchn(P, PQ);
      %  [idx(i), dist(i)] = dsearchn(P, P);
    %     figure,
    %     plot(P(:,1),P(:,2),'ko')
    %     hold on
    %     plot(PQ(:,1),PQ(:,2),'*g')
    %     hold on
    %     plot(P(idx(i),1),P(idx(i),2),'*r')
    %     legend('Data Points','Query Points','Nearest Points','Location','sw')

    end
    dist_all(j,:) = dist;
end
dist_all_2 = dist_all.*(wavelengths(k)./(2*pi))-2.*r;
%%
%mean(dist)
mean(dist_all,2).*(wavelengths(k)./(2*pi))
mean(dist_all(:)).*(wavelengths(k)./(2*pi))-2*r
mean(dist_all(:)).*(wavelengths(k)./(2*pi))./(2*r)

% parfor i = 1:lcords
%     
%     d(i) = dsearchn(lcords(i), lcords);
%     
%     
%     
%     
%     
% %    for j = 1:lcords     
% %        
% %        d = 
% % %         d(i, j) = get_distance(cords(j), cords(i));
% % %         if d(i,j) == 0
% % %             d(i,j) = 10000;
% % %         end
% %         
% %    end
% %    d_min(i) = min(d(i,:)
% end



function d = get_distance(c1, c2)
    x2 = (c1(1)-c2(1)).^2;
    y2 = (c1(2)-c2(2)).^2;
    z2 = (c1(3)-c2(3)).^2;
    d = sqrt(x2+y2+z2);
end


























