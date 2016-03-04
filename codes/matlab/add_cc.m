close all; 
clear all; 
clc;
if ~exist('s','var')
    prj_path='';
    proj_path;
    s = vml('cavriglia','2015_07_30',[],true);
    n_data={};
end
load('calc\data_with_sat_fact.mat', 'data');
all_j = cellfun(@(d) uint16(d.j), data);
all_days = cellfun(@(d) d.day, data,'UniformOutput', false);
last_batch = 1;
R = 60;
cl=1:300; 
rr=(1:300)';
mid_im_x=s.mfi.sz(1)/2;
mid_im_y=s.mfi.sz(2)/2;

t_l_sm = s.mfi.sm & repmat(cl<mid_im_y,s.mfi.sz(2),1) & repmat(rr<mid_im_x,1,s.mfi.sz(1));
t_r_sm = s.mfi.sm & repmat(cl>=mid_im_y,s.mfi.sz(2),1) & repmat(rr<mid_im_x,1,s.mfi.sz(1));
b_l_sm = s.mfi.sm & repmat(cl<mid_im_y,s.mfi.sz(2),1) & repmat(rr>=mid_im_x,1,s.mfi.sz(1));
b_r_sm = s.mfi.sm & repmat(cl>=mid_im_y,s.mfi.sz(2),1) & repmat(rr>=mid_im_x,1,s.mfi.sz(1));

t_l_area = sum(sum(t_l_sm));
t_r_area = sum(sum(t_r_sm));
b_l_area = sum(sum(b_l_sm));
b_r_area = sum(sum(b_r_sm));

while last_batch<30
    try
        if last_batch<11
            load(['D:\abb_data\cc_data',num2str(last_batch), '.mat'], 'ccp');
        else
            load(['D:\abb_data\cc_data',num2str(last_batch), '.mat'], 'cc');
            ccp=cc;
            clear cc;
        end
        disp(['processing cc_data',num2str(last_batch), '.mat ... ']);
        for i=1:length(ccp)    
            ccd = ccp{i};
            j_ind = all_j==ccd.j;
            ind_in_j = find(cellfun(@(dy) strcmp(dy,ccd.day), all_days(j_ind)));
            if isempty(ind_in_j), disp(['frame ', num2str(ccd.j),' on',ccd.day,' not found!']); continue; end
            ind = find(j_ind,ind_in_j,'first');
            d = data{ind(end)};
            disp(ind(end));
            clouds=false(300);
            clouds(ccd.ind)=true;
%             img_file = [img_save_path d.day '__' num2str(d.j) '.jpeg'];
%             I = s.get_image(img_file);
%             figure(2); imshow(I);
            sun_pos = s.sunpos_im(d.time).*(s.mfi.sz./s.oi.sz)';
            sx=floor(sun_pos(1)); sy=floor(sun_pos(2));
            f=@(xx,yy) (xx-sx).^2+(yy-sy).^2 <=R^2; 
            C=bsxfun(f,rr,cl);
            s_t_l = clouds & C & repmat(cl<sy,size(C,2),1) & repmat(rr<sx,1,size(C,1));
            s_t_r = clouds & C & repmat(cl>=sy,size(C,2),1) & repmat(rr<sx,1,size(C,1));
            s_b_l = clouds & C & repmat(cl<sy,size(C,2),1) & repmat(rr>=sx,1,size(C,1));
            s_b_r = clouds & C & repmat(cl>=sy,size(C,2),1) & repmat(rr>=sx,1,size(C,1));
            
            s_t_l_sm = s.mfi.sm & C & repmat(cl<sy,size(C,2),1) & repmat(rr<sx,1,size(C,1));
            s_t_r_sm = s.mfi.sm & C & repmat(cl>=sy,size(C,2),1) & repmat(rr<sx,1,size(C,1));
            s_b_l_sm = s.mfi.sm & C & repmat(cl<sy,size(C,2),1) & repmat(rr>=sx,1,size(C,1));
            s_b_r_sm = s.mfi.sm & C & repmat(cl>=sy,size(C,2),1) & repmat(rr>=sx,1,size(C,1));

            t_l = clouds & repmat(cl<mid_im_y,size(C,2),1) & repmat(rr<mid_im_x,1,size(C,1));
            t_r = clouds & repmat(cl>=mid_im_y,size(C,2),1) & repmat(rr<mid_im_x,1,size(C,1));
            b_l = clouds & repmat(cl<mid_im_y,size(C,2),1) & repmat(rr>=mid_im_x,1,size(C,1));
            b_r = clouds & repmat(cl>=mid_im_y,size(C,2),1) & repmat(rr>=mid_im_x,1,size(C,1));
            
            d.cc_fact = 100*[sum(sum(s_t_l))/sum(sum(s_t_l_sm)),sum(sum(s_t_r))/sum(sum(s_t_r_sm)),sum(sum(s_b_l))/sum(sum(s_b_l_sm)),sum(sum(s_b_r))/sum(sum(s_b_r_sm)),sum(sum(t_l))/t_l_area...
                sum(sum(t_r))/t_r_area,sum(sum(b_l))/b_l_area,sum(sum(b_r))/b_r_area];
            if isempty(ccd.ind)
                if d.clouds>90
                    d.cc_fact = repmat(d.clouds,1,8);
                    n_data{end+1} = d;
                end
            else
                n_data{end+1} = d;
            end
        end
    catch ME
        disp('AN ERROR!');
        disp(getReport(ME));
    end
    last_batch = last_batch + 1;
end
data=n_data;
save('calc\clean_data_with_8cc_nan_corrected.mat', 'data');
