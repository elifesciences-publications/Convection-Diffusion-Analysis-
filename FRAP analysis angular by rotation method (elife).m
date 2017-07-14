%##########################################################################
%###########     Convection Diffusion Analysis for Glymphatic  ############
%###########          edited by Byung-Ju Jin (2017)      ##################
%#############  University of California, San Francisco  ##################


clear all; 
cputime0 = cputime;

CROP = 1; % 0 - no crop, 1 - crop 
CONVEC = 0.; % artificial convection 10pixel * 0.4 um = 4um/s 
dn = 10; % add dn pixels between each pixel 

n_start = 1; % first file 
n_step = 1; % file interval 
n_end = 74; % last file 
n_file = 74; % total nubmer of file  

i_bleaching = 6; 

FPS = 0.77; % camera frame per second   
dtime = 1/FPS;


fname_head = '30%_200'; % file name
fname = [fname_head,'01.tif']; % initial file name 
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



xstep = 1.;  
ystep = 1.;  
xstep_n = 10*dn; 
ystep_n = 10*dn; 

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

[cx cy] = ginput(1);
[bx by] = ginput(1); 
radius = sqrt((cx-bx)^2. + (cy-by)^2.); 
t = linspace(0,2*pi,100);
x = radius*cos(t)+cx; 
y = radius*sin(t)+cy; 


N_ANGLE = 8; 
dANGLE = 2*pi/N_ANGLE; 
image_0 = zeros(n_row_new, n_col_new); 
mask_im = zeros(n_row_new, n_col_new,N_ANGLE); 
mask_im_cir = zeros(n_row_new, n_col_new); 
mask_im_sum = zeros(n_row_new,n_col_new,4) ; 

% determin the distance and angle 
for j=1:n_row_new
    for i=1:n_col_new
        theta = atan2((j-cy),(i-cx)); 
        r_dis = sqrt((j-cy)^2. + (i-cx)^2.);  
        
        if (r_dis < radius )
            for k_ANGLE=1:N_ANGLE
                ANGLE_min = -pi+dANGLE*(k_ANGLE-1); 
                ANGLE_max = ANGLE_min + dANGLE; 
                
                if ((theta >= ANGLE_min) && (theta <= ANGLE_max))
                    mask_im(j,i,k_ANGLE) = 1.; 
                end 
            end 
        end 
    end
end 


for j=1:n_row_new
    for i=1:n_col_new
        if ((mask_im(j,i,2) == 1) || (mask_im(j,i,3) == 1))
            mask_im_sum(j,i,1) = 1; 
        elseif ((mask_im(j,i,4) == 1) || (mask_im(j,i,5) == 1)) 
            mask_im_sum(j,i,2) = 1; 
        elseif ((mask_im(j,i,6) == 1) || (mask_im(j,i,7) == 1)) 
            mask_im_sum(j,i,3) = 1; 
        elseif ((mask_im(j,i,8) == 1) || (mask_im(j,i,1) == 1))             
            mask_im_sum(j,i,4) = 1; 
        else
            ;
        end
    end
end

%######## multiple file analysis  ####
kfile = 0; 

CONVEC = 1; % convective movement 
for ifile = n_start:n_step:n_end   

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

    for k_ANGLE=1:4
        imsi_mask = mask_im_sum(:,:,k_ANGLE); 
        mean_intensity_trace(kfile,k_ANGLE) = mean(mean(crop_im_new.*imsi_mask(:,:)));
    end 

end

for k_ANGLE=1:4
    mean_intensity_trace_norm(:,k_ANGLE) = mean_intensity_trace(:,k_ANGLE)/mean_intensity_trace(1,k_ANGLE)
end
    

    
    
    
    
    
    