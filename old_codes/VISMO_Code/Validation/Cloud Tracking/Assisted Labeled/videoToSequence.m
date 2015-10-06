filename = 'Validation_Sequence_231.avi';

mov = VideoReader(filename);
iterator = 1;

for i = 7:4:42
    ff(iterator) = {read(mov,i)};
    iterator = iterator +1;
end
row_c = 3;
col_c = 3;
result_im = zeros(size(cell2mat(ff(1))).*[row_c,col_c,1]);

for i = 1:length(ff)
    frame = cell2mat(ff(i));
    r_c = size(frame,1);
    c_c = size(frame,2);
    r = double(idivide((i-1),uint8(col_c))+1);
    c = mod(i-1,col_c)+1;
    result_im(((r-1)*r_c)+1:r*r_c,((c-1)*c_c)+1:c*c_c,:) = uint8(frame);    
end

imshow(uint8(result_im));
imwrite(uint8(result_im),'sequence_266.jpg');