function lut = getLUT(sph_im,opt,ocam_model)
imsize = size(sph_im);

% extracting the usefull parts of the image
margin_f = 0;
image_width = imsize(1)-2*margin_f;

% the actual center of the real image.
deg2rad = pi/180;
cent = ceil(imsize(1:2)/2);
fov = 185*deg2rad;
k = 1.2;
f_p = image_width/(4*sin(fov/4));
out_cent = ceil(cent*0.3);
%out_cent = cent;
%image_width
%out_cent = [ round(f_p * tan( asin( sin(atan((image_width/2)/f_p))*k ) )), round(f_p * tan( asin( sin(atan((image_width/2)/f_p))*k ) ))];
%str_cent = [ round(f_p * tan( asin( sin(atan(((image_width-2*margin_f)/2)/f_p))*k ) )), round(f_p * tan( asin( sin(atan(((image_width-2*margin_f)/2)/f_p))*k ) ))];
sph_im = double(sph_im);

out_size = numel(-(out_cent(1)-1):(out_cent(1)-1));

lut = struct(...
    'valid', {},...
    'pnnx', {},...
    'pnny', {});
tic
% c = ocam_model.c;
% d = ocam_model.d;
% e = ocam_model.e;
% pol_params = ocam_model.ss;
%
if(opt == 1)
    %s = out_cent(1)/str_cent(1);
    
    % for a square image this should be the case...
    for ij=-(out_cent(1)-1):(out_cent(1)-1)
        %cf_ij = (out_cent(1)-abs(ij));
        %cf_ij = 10.2;
        %m_ij(ij+out_cent(1)) = double(ij);
        %mm_ij(ij+out_cent(1)) = cf_ij;
        %mm_ij(ij+out_cent(1)) = atan(double(ij)/cf_ij);
        %mmm_ij(ij+out_cent(1)) = double(sin( atan(ij/cf_ij) )) / s ;
        %mmmm_ij(ij+out_cent(1)) = asin( sin( atan(ij/cf_ij) ) / s );
        %map_ij(ij+out_cent(1)) =  cf_ij*tan( asin( (sin( atan(ij/cf_ij) ) / s) ) );
        map_ij(ij+out_cent(1)) = ij;
    end
    %plot(-(out_cent(1)-1):(out_cent(1)-1),map_ij)
    %hold on
    %plot(-(out_cent(1)-1):(out_cent(1)-1),m_ij,'g')
    %plot(-(out_cent(1)-1):(out_cent(1)-1),mm_ij,'k')
    
    for i=-(out_cent(1)-1):(out_cent(1)-1)
        for j=-(out_cent(2)-1):(out_cent(2)-1)
            %lambda = sqrt(i^2+j^2)/f_p;
            %m = world2cam([i,j,-sqrt((out_cent(1)-1)^2-i^2-j^2)]',ocam_model);
            m = world2cam([i,j,-10]',ocam_model);
            x_f = (m(1)-823)*(823/823);
            y_f = (m(2)-823)*(823/823);
            %x_f = 2*i*sin(atan(lambda)/2)/lambda;
            %y_f = 2*j*sin(atan(lambda)/2)/lambda;
            
            %m = world2cam([i,j,1]',ocam_model);
            %x_f = m(1);
            %y_f = m(2);
            
            %[theta,R] = cart2pol(map_ij(i+out_cent(1)),map_ij(j+out_cent(2)));
            %s = pol_params(1) + pol_params(2)*r + pol_params(3)*r.^2 +pol_params(4)*r.^3 + pol_params(5)*r.^4;
            %r = f_p * tan( asin( sin(atan(R/f_p))/k ) );
            %[x_f,y_f] = pol2cart(theta,r);
            %x_f = c*i + d*j;
            %y_f = e*i + j;
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
                
                if( (max(pnnx(:)) <= image_width) && (max(pnny(:)) <= image_width) && (min(pnnx(:)) > 0) && (min(pnny(:)) > 0))
                    new_str = struct(...
                        'valid', 1,...
                        'pnnx', pnnx,...
                        'pnny', pnny);
                    lut(i+out_cent(1),j+out_cent(2)) = new_str;
                    %rect_im(i+out_cent(1),j+out_cent(2),:) = sum(sum(sph_im(pnnx(:),pnny(:),:),1),2)/(numel(pnnx(:))*numel(pnny(:)));
                else
                    pnnx = 0;
                    pnny = 0;
                    new_str = struct(...
                        'valid', 0,...
                        'pnnx', pnnx,...
                        'pnny', pnny);
                    lut(i+out_cent(1),j+out_cent(2)) = new_str;
                    %rect_im(i+out_cent(1),j+out_cent(2),:) = [0,0,0];
                end
                %rect_im(i+out_cent(1),j+out_cent(2),:) = sph_im(ceil(x_f+cent(1)),ceil(y_f+cent(2)),:);
            else
                pnnx = 0;
                pnny = 0;
                new_str = struct(...
                    'valid', 0,...
                    'pnnx', pnnx,...
                    'pnny', pnny);
                lut(i+out_cent(1),j+out_cent(2)) = new_str;
                %rect_im(i+out_cent(1),j+out_cent(2),:) = [0,0,0];
            end
        end
    end
elseif(opt == 2)
    margin_f = 17;
    image_width = imsize(1)-2*margin_f;
    out_cent = ceil(cent*1.4);
    f_p = image_width/(4*sin(fov/2));
    out_size = numel(-(out_cent(1)-1):(out_cent(1)-1));
    % for a square image this should be the case...
    for i=-(out_cent(1)-1):(out_cent(1)-1)
        for j=-(out_cent(2)-1):(out_cent(2)-1)
            lambda = sqrt(i^2+j^2)/f_p;
            x_f = 2*i*sin(atan(lambda)/2)/lambda;
            y_f = 2*j*sin(atan(lambda)/2)/lambda;
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
                
                if( (max(pnnx(:)) <= image_width) && (max(pnny(:)) <= image_width) && (min(pnnx(:)) > 0) && (min(pnny(:)) > 0))
                    new_str = struct(...
                        'valid', 1,...
                        'pnnx', pnnx,...
                        'pnny', pnny);
                    lut(i+out_cent(1),j+out_cent(2)) = new_str;
                    %rect_im(i+out_cent(1),j+out_cent(2),:) = sum(sum(sph_im(pnnx(:),pnny(:),:),1),2)/(numel(pnnx(:))*numel(pnny(:)));
                else
                    pnnx = 0;
                    pnny = 0;
                    new_str = struct(...
                        'valid', 0,...
                        'pnnx', pnnx,...
                        'pnny', pnny);
                    lut(i+out_cent(1),j+out_cent(2)) = new_str;
                    %rect_im(i+out_cent(1),j+out_cent(2),:) = [0,0,0];
                end
                %rect_im(i+out_cent(1),j+out_cent(2),:) = sph_im(ceil(x_f+cent(1)),ceil(y_f+cent(2)),:);
            else
                pnnx = 0;
                pnny = 0;
                new_str = struct(...
                    'valid', 0,...
                    'pnnx', pnnx,...
                    'pnny', pnny);
                lut(i+out_cent(1),j+out_cent(2)) = new_str;
                %rect_im(i+out_cent(1),j+out_cent(2),:) = [0,0,0];
            end
        end
    end
else
    [out_y,out_x] = meshgrid(-(out_cent(1)-1):(out_cent(1)-1),-(out_cent(2)-1):(out_cent(2)-1));
    lut = cell2mat(arrayfun(@(x,y) getCorrPt(sph_im,cent,out_size,f_p,x,y),out_x,out_y,'UniformOutput',false));
end
toc
end


function lut_pt = getCorrPt(sph_im,cent,out_size,f_p,i,j)

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
        lut_pt = struct(...
            'valid', 1,...
            'pnnx', pnnx,...
            'pnny', pnny);
    else
        lut_pt = struct(...
            'valid', 0,...
            'pnnx', {},...
            'pnny', {});
    end
    
    % no interpolation, near neighbour mapping...
    %rgb_val(1,1,:) = sph_im(ceil(x_f+cent(1)),ceil(y_f+cent(2)),:);
else
    lut_pt = struct(...
        'valid', 0,...
        'pnnx', {},...
        'pnny', {});
end

end


