%% This tool is prepared to be able to annotate on images using different 
% shapes and marker styles. Optimized for rectangle plots.
% B. Zeydan, 27. Feb. 2013

function im = objectAnnotation(im,type,ants,color,thickness)
if(size(im,3)<3)
    [im,im_map] = gray2ind(im);
    im = ind2rgb(im,im_map);
end

if nargin < 4
    color = 'y';
    thickness = 1;
elseif nargin < 5
    thickness = 1;
end

switch char(color)
    case 'r'
        c = [255,-255,-255];
    case 'g'
        c = [-255,255,-255];
    case 'b'
        c = [-255,-255,255];
    case 'y'
        c = [255,255,-255];
    case 'c'
        c = [-255,255,255];
    case 'm'
        c = [255,-255,255];
    case 'p'
        c = [255,-255,255];
    case 'k'
        c = [-255,-255,-255];
    case 'w'
        c = [255, 255, 255];
    otherwise
        inds = color < 255;
        color(inds) = (-255*ones(size(color(inds))))+color(inds);
        c = color;
end

switch type
    case 'Rectangle'
        bx_pts = cell(size(ants,1)*thickness,1);
        count = 1;
        for i = 1:size(ants,1)
            
            o_bbox = ants(i,:);
            for j = 0:thickness-1
                bbox = o_bbox + [j j -2*j -2*j];
                bx_pts{count} = [repmat(bbox(2),length(bbox(1):bbox(1)+bbox(3)),1),[bbox(1):bbox(1)+bbox(3)]';...
                    repmat(bbox(2)+bbox(4),length(bbox(1):bbox(1)+bbox(3)),1),[bbox(1):bbox(1)+bbox(3)]';...
                    [bbox(2):bbox(2)+bbox(4)]',repmat(bbox(1),length(bbox(2):bbox(2)+bbox(4)),1);...
                    [bbox(2):bbox(2)+bbox(4)]',repmat(bbox(1)+bbox(3),length(bbox(2):bbox(2)+bbox(4)),1)];
                count = count + 1;
            end
            
        end
        
        bx_pts = cell2mat(bx_pts(:));
        valid_x = (bx_pts(:,1) <= size(im,1)) & (bx_pts(:,1) >= 1);
        valid_y = (bx_pts(:,2) <= size(im,2)) & (bx_pts(:,2) >= 1);
        valid_pts = valid_x & valid_y;
        
        bx_pts = bx_pts(valid_pts,:);
        ibx_pts = sub2ind(size(im),bx_pts(:,1),bx_pts(:,2));
        
        mask = zeros(size(im,1),size(im,2));
        mask(ibx_pts) = 1;
        mask3 = cat(3,c(1)*mask,c(2)*mask,c(3)*mask);
        
        im = im + double(mask3);
    case 'CVRectangles'
        for i = 1:size(ants,1)
            bbox = ants(i,:);
            if(bbox(3)<1 || bbox(4)<1)
                continue;
            end
            im = cv.rectangle(im,bbox,'Color',c);
        end
    case 'CVArrow'
        for i = 1:size(ants,1)
            im = cv.line(im,[ants(i,1),ants(i,2)]-1,[ants(i,3),ants(i,4)]-1,'Color',c,'LineType',4,'Thickness',1.8);
        end
    case 'Text'
        for i = 1:size(ants,1)
            im = cv.putText(im,num2str(ants(i,1)),[ants(i,2),ants(i,3)],'Color',c , 'Thickness', thickness);
        end
    otherwise
        disp('Annotation type not recognized. Try Rectangle or Label...');
end
end