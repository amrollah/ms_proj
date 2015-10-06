%% Circular Fish Eye Defishing Code with 9-Point interpolation.
%
% standard Defishing implementation for could tracking:
% defish_im: outpu, defished resulting image
% fish_im: input, circular fish eye image
%
% B. Zeydan, 05. Mar. 2013

function [ rect_im ] = sph2rect( sph_im , lut)

if(nargin < 2)
    imsize = size(sph_im);
    
    % extracting the usefull parts of the image
    margin = 17;
    image_width = imsize(1)-2*margin;
    
    % the actual center of the real image.
    cent = ceil(imsize(1:2)/2);
    out_cent = ceil(cent*1.4);
    
    deg2rad = 0.0174532925;
    fov = 185*deg2rad;
    f_p = image_width/(4*sin(fov/2));
    sph_im = double(sph_im);
    
    out_size = numel(-(out_cent(1)-1):(out_cent(1)-1));
    
    [out_y,out_x] = meshgrid(-(out_cent(1)-1):(out_cent(1)-1),-(out_cent(2)-1):(out_cent(2)-1));
    rect_im = cell2mat(arrayfun(@(x,y) getRGB(sph_im,cent,out_size,f_p,x,y),out_x,out_y,'UniformOutput',false));  
else
    %out_y,out_x] = meshgrid(1:size(lut,1),1:size(lut,2));
    %rect_im = cell2mat(arrayfun(@(x,y) lutRGB(sph_im,lut,x,y),out_x,out_y,'UniformOutput',false));
    rect_im = zeros(size(lut,1),size(lut,2),3);
    for k = 1:size(lut,1)
        for l = 1:size(lut,2)
            if((lut(k,l).valid == 0))
                rect_im(k,l,:) = reshape([0,0,0],[1,1,3]);
            else
                pnnx = lut(k,l).pnnx;
                pnny = lut(k,l).pnny;
                rect_im(k,l,:) = sum(sum(sph_im(pnnx(:),pnny(:),:),1),2)/(numel(pnnx(:))*numel(pnny(:)));
            end
        end
    end
end
rect_im = uint8(rect_im);
end

function lut = getLUT(sph_im)
imsize = size(sph_im);

% extracting the usefull parts of the image
margin = 17;
image_width = imsize(1)-2*margin;

% the actual center of the real image.
cent = ceil(imsize(1:2)/2);
out_cent = ceil(cent*1.4);

deg2rad = 0.0174532925;
fov = 185*deg2rad;
f_p = image_width/(4*sin(fov/2));
sph_im = double(sph_im);

out_size = numel(-(out_cent(1)-1):(out_cent(1)-1));

lut = cell(length(-(out_cent(1)-1):(out_cent(1)-1)),length(-(out_cent(1)-1):(out_cent(1)-1)));

for i=-(out_cent(1)-1):(out_cent(1)-1)
    for j=-(out_cent(2)-1):(out_cent(2)-1)
        x_f = 2*i*sin(atan(sqrt(i^2+j^2)/f_p))/(sqrt(i^2+j^2)/f_p);
        y_f = 2*j*sin(atan(sqrt(i^2+j^2)/f_p))/(sqrt(i^2+j^2)/f_p);
        if((x_f+cent(1)) > 0 && (y_f+cent(2)) > 0)
            % floor(x_f+cent(1)):0.5:ceil(x_f+cent(1)),floor(y_f+cent(2)):0.5:ceil(y_f+cent(2)) grid of 9-near neighboor pixels
            % we need to find the closest and give its value to the result
            [nnx,nny] = meshgrid(floor(x_f):0.5:ceil(x_f),floor(y_f):0.5:ceil(y_f));
            comp_g = [nnx(:),nny(:)];
            comp_p = repmat([x_f,y_f],numel(comp_g(:,1)),1);
            comp_r = comp_g-comp_p;
            [~,ires_p] = min(sqrt(comp_r(:,1).^2 + comp_r(:,2).^2));
            res_p = comp_g(ires_p(1),:);
            %[pnnx,pnny] = meshgrid(floor(res_p(1)+cent(1)):ceil(res_p(1)+cent(1)),floor(res_p(2)+cent(2)):ceil(res_p(2)+cent(2)));
            pnnx = floor(res_p(1)+cent(1)):ceil(res_p(1)+cent(1));
            pnny = floor(res_p(2)+cent(2)):ceil(res_p(2)+cent(2));
            
            if( (max(pnnx(:)) <= out_size) && (max(pnny(:)) <= out_size) )
                lut(i+out_cent(1),j+out_cent(2)) = [pnnx(:),pnny(:)];
                %rect_im(i+out_cent(1),j+out_cent(2),:) = sum(sum(sph_im(pnnx(:),pnny(:),:),1),2)/(numel(pnnx(:))*numel(pnny(:)));
            else
                lut(i+out_cent(1),j+out_cent(2)) = 0;
                %rect_im(i+out_cent(1),j+out_cent(2),:) = [0,0,0];
            end
            %rect_im(i+out_cent(1),j+out_cent(2),:) = sph_im(ceil(x_f+cent(1)),ceil(y_f+cent(2)),:);
        else
            lut(i+out_cent(1),j+out_cent(2)) = 0;
            %rect_im(i+out_cent(1),j+out_cent(2),:) = [0,0,0];
        end
    end
end
end

function rgb_val = getRGB(sph_im,cent,out_size,f_p,i,j)
rgb_val = zeros(1,1,3);
lambda = sqrt(i^2+j^2)/f_p;
x_f = 2*i*sin(atan(lambda/2))/lambda;
y_f = 2*j*sin(atan(lambda/2))/lambda;

if((i == 0) && (j == 0)) % to handle the center point...
    x_f = 0;
    y_f = 0;
end

if(((x_f+cent(1)) > 0 && (y_f+cent(2)) > 0))  % the last is the center condition to be captured...
    
    % the part to do the interpolation. If the interpolation is not
    % required, than skip this part as it slows down the method...
    
    % floor(x_f+cent(1)):0.5:ceil(x_f+cent(1)),floor(y_f+cent(2)):0.5:ceil(y_f+cent(2)) grid of 9-near neighboor pixels
    % we need to find the closest and give its value to the result
    [nnx,nny] = meshgrid(floor(x_f):0.5:ceil(x_f),floor(y_f):0.5:ceil(y_f));
    comp_g = [nnx(:),nny(:)];
    comp_p = repmat([x_f,y_f],numel(comp_g(:,1)),1);
    comp_r = comp_g-comp_p;
    [~,ires_p] = min(sqrt(comp_r(:,1).^2 + comp_r(:,2).^2));
    res_p = comp_g(ires_p(1),:);
    pnnx = floor(res_p(1)+cent(1)):ceil(res_p(1)+cent(1));
    pnny = floor(res_p(2)+cent(2)):ceil(res_p(2)+cent(2));
    
    if( (max(pnnx(:)) <= out_size) && (max(pnny(:)) <= out_size) )
        rgb_val(1,1,:) = sum(sum(sph_im(pnnx(:),pnny(:),:),1),2)/(numel(pnnx(:))*numel(pnny(:)));
    else
        rgb_val(1,1,:) = [0,0,0];
    end
    
    % no interpolation, near neighbour mapping...
    %rgb_val(1,1,:) = sph_im(ceil(x_f+cent(1)),ceil(y_f+cent(2)),:);
else
    rgb_val(1,1,:) = [0,0,0];
end

end

function rgb_val = lutRGB(sph_im,lut,k,l)

if((lut(k,l).valid == 0))
    rgb_val = reshape([0,0,0],[1,1,3]);
else
    pnnx = lut(k,l).pnnx;
    pnny = lut(k,l).pnny;
    rgb_val = sum(sum(sph_im(pnnx(:),pnny(:),:),1),2)/(numel(pnnx(:))*numel(pnny(:)));
end
end
