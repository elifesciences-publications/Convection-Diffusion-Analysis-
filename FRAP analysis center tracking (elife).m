%##########################################################################
%###########     Convection Diffusion Analysis for Glymphatic  ############
%###########          edited by Byung-Ju Jin (2017)      ##################
%#############  University of California, San Francisco  ##################


clear all; 
cputime0 = cputime;

CROP = 1; % 0 - no crop, 1 - crop 
CONVEC = 0.0; % artificial convection 10pixel * 0.4 um = 4um/s 
dn = 5; % add dn pixels between each pixel 

% file parameter set 
n_start = 1; % first file (initial bleaching)
n_step = 1; % file interval  
n_end = 74; % last file 
n_file = n_end - n_start +1; % total nubmer of file  

i_bleaching = 7; % the beginning file of the bleaching 
% FPS = 0.77; % camera frame per second   
% dtime = 1/FPS;

pixeldimension = 0.4; % um 
pixeldimension_new = pixeldimension/dn; 

% sweeping parameter set 
xstep = 1.; % sweeping step (pixel)
ystep = 1.;  
xstep_n = 10*dn; % sweeping range (pixel) 
ystep_n = 10*dn; 

% file reading 
fname_head = '30%_200'; % file name
fname = [fname_head,'06.tif']; % initial file name 
image = double(imread(fname,'tif'));
imshow(image/max(max(image))); 


if (CROP == 1)
    sprintf('crop the ROI ----->') 
    [a b] = ginput(2); 
    x_cr=round(a);
    y_cr=round(b);
    crop_im = image(y_cr(1):y_cr(2),x_cr(1):x_cr(2));
    [n_row n_col] = size(crop_im); 
else 
    crop_im = image; 
    [n_row n_col] = size(crop_im); 
end

% add pixels to increase the resolution 
for j=1:(n_row-1)
    for i=1:(n_col-1)
        
        for m=1:(dn+1)
            for n=1:(dn+1) 
                crop_im_new((j+m-1+dn*(j-1)),(i+n-1+dn*(i-1))) = crop_im(j,i); 
            end
        end
    end
end


crop_im_new((dn*(n_row-1)+n_row), (dn*(n_col-1)+n_col))= crop_im(n_row,n_col); 
crop_im_new_save = crop_im_new; 

n_row_new = dn*(n_row-1)+n_row; 
n_col_new = dn*(n_col-1)+n_col; 


imshow(crop_im_new/max(max(crop_im_new))); 
% figure, imshow(crop_im/max(max(crop_im))); 


mask_im = zeros(n_row_new, n_col_new); % Image mask of bleaching area  

sprintf('point the center point and the boundary point----->')
[cx cy] = ginput(1);
[bx by] = ginput(1); 
radius = sqrt((cx-bx)^2. + (cy-by)^2.); 
t = linspace(0,2*pi,1000);
x = radius*cos(t)+cx; 
y = radius*sin(t)+cy; 

hold on 
plot(x,y,'b','LineWidth',2); 


k = 0; 
for j=1:n_row_new
    for i=1:n_col_new
        distance = sqrt((cx-i)^2. + (cy-j)^2.); 
        if (distance < radius) 
            k = k +1; 
            in_x(k) = i;
            in_y(k) = j; 
            mask_im(j,i) = 1.; 
        else
            ;
        end 
    end
end
kn = k; % number of points in bleaching area 
% sweep in x, y direction and calculate the mean intensity 

kmask = 0; 
for n=-xstep_n:xstep:xstep_n
    for m=-ystep_n:ystep:ystep_n
        mask_im = zeros(n_row_new, n_col_new); 
        kmask = kmask + 1; 
        for k=1:kn
            mask_im(in_y(k)+m,in_x(k)+n) = 1.; 
        end
        mean_intensity(ystep_n+m+1,xstep_n+n+1) = mean(mean(crop_im_new.*mask_im)); 
    end
end


minimum_intensity = min(min(mean_intensity)); 
[dy, dx] = find(mean_intensity == minimum_intensity); 
if (size(dx) == 1) 
    ycenter = cy - (ystep_n - dy*ystep); 
    xcenter = cx - (xstep_n - dx*xstep); 
else 
    dx_mean = mean(dx); 
    dy_mean = mean(dy); 
    ycenter = cy - (ystep_n - dy_mean*ystep); 
    xcenter = cx - (xstep_n - dx_mean*xstep); 
end 

% confirm 
t = linspace(0,2*pi,1000);
x = radius*cos(t)+xcenter; 
y = radius*sin(t)+ycenter; 

hold on 
plot(x,y,'r','LineWidth',2); 


% Multiple image analysis 
%######## multiple file analysis  ####
kfile = 0; 

for ifile = n_start:n_step:n_end   
% for ifile = n_start:n_step:11   
    fname = [fname_head, sprintf('%02d.tif',ifile)]; % 02d -> 00, 01, 02, 03... 
                                                     % 03d -> 001, 002, 003
    kfile = kfile+1; 

    k_image = double(imread(fname,'TIF'));

    if (CROP == 1)
        if (ifile > i_bleaching)
            movement = round(CONVEC*(ifile-i_bleaching))
            crop_im = k_image(y_cr(1):y_cr(2),(x_cr(1)-movement):(x_cr(2)-movement));
            ifile 
        else
            crop_im = k_image(y_cr(1):y_cr(2),x_cr(1):x_cr(2));
        end 
    else 
        crop_im = k_image; 
    end
    
    % add pixels to increase the resolution 

    for j=1:(n_row-1)
        for i=1:(n_col-1)

            for m=1:(dn+1)
                for n=1:(dn+1) 
                    crop_im_new((j+m-1+dn*(j-1)),(i+n-1+dn*(i-1))) = crop_im(j,i); 
                end
            end
        end
    end

    

    crop_im_new((dn*(n_row-1)+n_row), (dn*(n_col-1)+n_col))= crop_im(n_row,n_col); 



% sweep in x, y direction and calculate the mean intensity 

    kmask = 0; 
    for n=-xstep_n:xstep:xstep_n
        for m=-ystep_n:ystep:ystep_n
            mask_im = zeros(n_row_new,n_col_new); 
            kmask = kmask + 1; 
            for k=1:kn
                mask_im(in_y(k)+m,in_x(k)+n) = 1.; % move the mask m by n distance 
            end
            mean_intensity_trace(ystep_n+m+1,xstep_n+n+1,kfile) = mean(mean(crop_im_new.*mask_im)); 
            mean_intensity_trace_norm(ystep_n+m+1,xstep_n+n+1,kfile) ...
                = mean_intensity_trace(ystep_n+m+1,xstep_n+n+1,kfile)/...
                mean_intensity_trace(ystep_n+m+1,xstep_n+n+1,1); 
        end
    end

end
crop_im_final = crop_im_new; 
imshow(crop_im_final/max(max(crop_im_final))); 

for i = 1:kfile
   minimum_value = min(min(mean_intensity_trace_norm(:,:,i)));
   [dy, dx] = find(mean_intensity_trace_norm(:,:,i)==minimum_value); 
   if (size(dx) == 1)
       xcenter_final_trace(i) = cx - (xstep_n - dx*xstep);
       ycenter_final_trace(i) = cy - (ystep_n - dy*ystep); 
   else 
       xcenter_final_trace(i) = cx - (xstep_n - mean(dx)*xstep);
       ycenter_final_trace(i) = cy - (ystep_n - mean(dy)*ystep); 
   end 

end 

   
% distance from the initial cente
for i = 1:kfile
       center_distance(i) = pixeldimension_new*sqrt((xcenter_final_trace(i) - xcenter_final_trace(i_bleaching))^2....
           +(ycenter_final_trace(i) - ycenter_final_trace(i_bleaching))^2.);  
end 


   
sprintf('bleaching area radius (um) ----->'); 
radius*pixeldimension
sprintf('maximum displacement from the center (pixels)---->'); 
max(center_distance(i_bleaching:(kfile-1)))
sprintf('mean displacement from the center (pixels)---->'); 
mean(center_distance(i_bleaching:(kfile-1)))
% i_bleaching: the beginning file of the bleaching 

% to plot recovery at three different points 
for i=1:kfile
    plotpoints(i,1) = mean_intensity_trace_norm(ystep_n, xstep_n,i); 
    plotpoints(i,2) = mean_intensity_trace_norm(ystep_n, xstep_n+round(radius),i); 
    plotpoints(i,3) = mean_intensity_trace_norm(ystep_n, xstep_n+2*round(radius),i); 
end 



% Image processing for journal figures
imshow(crop_im_new_save/max(max(crop_im_new_save))); 
hold on 


% initial circle by manual 
t = linspace(0,2*pi,1000);
x = radius*cos(t)+cx; 
y = radius*sin(t)+cy; 

plot(x,y,'--r','LineWidth',2); 

            
hold on
plot(xcenter_final_trace(i_bleaching:(kfile-1)),...
    ycenter_final_trace(i_bleaching:(kfile-1)),...
    '-ys','LineWidth',2,...
                'MarkerEdgeColor','g',...
                'MarkerFaceColor','g',...
                'MarkerSize',3);

