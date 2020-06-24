clc;

r = 44; 
[cords]  = get_cords(sphere_coeffs_perp, sphere_coeffs_par, r);
%cords = cords(:,:,1); %+(20*r+2*r)/2;
I = zeros(20*r+2*r,20*r+2*r);
I = draw_spheres(I, cords, r);

figure, imshow(I(:,:,1))



cors = get_correlations(I);



%%
clc
figure,
dk = 1./size(I,1);
fontSize = 24;
kx = ((1:1:size(cors,1))-size(cors,1)./2).*dk;
ky = ((1:1:size(cors,2))-size(cors,2)./2).*dk;
imshow(mean(cors,3), 'XData',kx , 'YData', ky);
ylim([-0.08,0.08])
xlim([-0.08,0.08])
%title('Cosine Function, f=(4,3)c/ph', 'FontSize', fontSize)
xlabel('k_x (nm^{-1})', 'FontSize', fontSize)
ylabel('k_y (nm^{-1})', 'FontSize', fontSize)
% Put tick marks at specified locations:
axis on
cb = colorbar;
ylabel(cb, 'Correlation coefficient')
xticks([-0.08:0.04:0.08])
yticks([-0.08:0.04:0.08])
h = gca;
set(h, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.85, 0.85]);
pbaspect([1 1 1])

set(gca,'FontSize',34)
title('30% Fill fraction') 
%%
% % Put lines across the entire image, in the overlay, where the main axes should be.
% line(xlim, [0,0], 'LineWidth', 2, 'Color', 'r');
% line([0,0], ylim, 'LineWidth', 2, 'Color', 'r');
len = size(cors,1);
line_profile = mean(cors(round(len/2),:,:),3);
figure, plot(kx,line_profile, 'LineWidth',8)
xlabel('k_x (nm^{-1})')
ylabel('Correlation coefficient')
set(gca,'FontSize',34)
title('30% Fill fraction (88nm Si)')
pbaspect([1 1 1])
ylim([-0.2,1.05])



%%
%figure, imshow(I(:,:,2))

[X,Y,Z] = sphere(280);
r = 44;
instance = 4;
figure, 
%h = light;
for idx = 1:size(cords,2)
    h = surf(X*r-cords(1,idx, instance),Y*r-cords(2,idx, instance),Z*r, 'EdgeColor', 'none', 'FaceColor',[0, 0.4470, 0.7410]);
    hold on 
    %lighting gouraud
    lightangle(-30,70)
    h.DiffuseStrength = 1;
    
    %camlight('right')
end

pbaspect([1 1 1])






% figure, imshow(mean(cors,3))


% figure, 
% imagesc(-12:1:12,-1:1:1,(mean(cat(3,cperp,cpar),3)))
% ylabel('Lag (Dimension)')
% xlabel('Lag (Particle)')
% colormap jet
% h = colorbar;
% ylabel(h,'Average Correlation Coefficient');
% yticks([-1 0 1])
% pbaspect([1 1 1])
% set(gca,'FontSize',34)
% title('20% Fill fraction') 


%%

function [all] = get_cords(sphere_coeffs_perp, sphere_coeffs_par, r)

for idx1 = 1:size(sphere_coeffs_perp,2)
    x_perp = [];
    y_perp = [];
    x_par = [];
    y_par = [];
    dummy_perp = sphere_coeffs_perp{1,idx1};
    dummy_par = sphere_coeffs_par{1,idx1};    
    ka = dummy_perp.ka;
    for idx2 = 1:size(dummy_perp, 2)
        x_perp = [x_perp,dummy_perp(idx2).loc_x];
        x_par = [x_par,-dummy_par(idx2).loc_y];
        y_perp = [y_perp,dummy_perp(idx2).loc_y];
        y_par = [y_par,dummy_par(idx2).loc_x];       
    end

    perp(:,:,idx1) = [x_perp;y_perp];
    par(:,:,idx1) = [x_par;y_par];
    allx = [x_perp,x_par];
    ally = [y_perp, y_par];
    all(:,:,idx1) = ([x_perp ; y_perp].*r./ka); %+11.*r;
%     call(:,:,idx1) = xcorr2(all(:,:,idx1));
%     cperp(:,:,idx1)  = normxcorr2(perp(:,:,idx1), perp(:,:,idx1));
%     cpar(:,:,idx1)  = normxcorr2(par(:,:,idx1), par(:,:,idx1) );

end
end

 
function I2 = draw_spheres(I, cords, r)
for idx1 = 1:size(cords, 3)
    dummy = I;
    for idx2 = 1:size(cords, 2)
        circle = make_circle(I, cords(:,idx2,idx1), r);
        dummy = dummy+circle;
    end 
    I2(:,:,idx1) = dummy;
end




end


function circle = make_circle(I, cord, r)

imageSizeX = size(I,1);
imageSizeY = size(I,2);
[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);

columnsInImage = columnsInImage-imageSizeX./2;
rowsInImage = rowsInImage-imageSizeY./2;

% Next create the circle in the image.
centerX = cord(1);
centerY = cord(2);
radius = r;
circle= (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;

end


function cors = get_correlations(I)

for idx = 1:1:size(I,3)
    cors(:,:,idx) = normxcorr2(I(:,:,idx),I(:,:,idx));
end

end






% function image = make_circle(cords, r, image)
% 
% 
% 
% 
% 
% 
% % Create a logical image of a circle with specified
% % diameter, center, and image size.
% % First create the image.
% imageSizeX = 640;
% imageSizeY = 480;
% [columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% % Next create the circle in the image.
% centerX = 320;
% centerY = 240;
% radius = 100;
% circlePixels = (rowsInImage - centerY).^2 ...
%     + (columnsInImage - centerX).^2 <= radius.^2;
% % circlePixels is a 2D "logical" array.
% % Now, display it.
% image(circlePixels) ;
% colormap([0 0 0; 1 1 1]);
% title('Binary image of a circle');

