function [mycc]=getcolours(myc,figON)
% function [mycc]=getcolours(myc,figON)
% Required input:  myc  can be either a number between 1 and 40
%                       corresponding to the official ABB colours or a
%                       string out of the following set:
%                       {'red','grey1','grey2','grey3','grey4','grey5','whi
%                       te','black'}. If myc exceeds 40 the colour rgb
%                       value will correspond to 
%                       rem(myc,40) and thus the function will cycle
%                       through the color codes.
% Optional input: figON If input is 1 then a figure with patch showing all
%                       possible colors will be plotted.
% Output:         mycc  The color code vector with the corresponding rgb
%                       values.
mycc = NaN;
if nargin < 2, figON = 0; end

c = defclr;
if figON,
    plotclrpatch(c);
end;

if ~ischar(myc),
    lc = length(c.C);
    r = abs(rem(myc,lc));
    if r == 0, 
        r=lc;
    end;
    mycc = c.C{r};
else
    
    switch lower(myc)
    case {'red','abbred','abb'},
        mycc = c.red;
    case {'grey1','abbgrey'},
        mycc = c.grey{1};
    case 'grey2',
        mycc = c.grey{2};
    case 'grey3',
        mycc = c.grey{3};
    case 'grey4',
        mycc = c.grey{4};
    case 'grey5',
        mycc = c.grey{5};
    case 'white' % 'white'
        mycc = c.white;
    case 'black' % 'black'
        mycc = c.black;
    otherwise % 'ABB red'
        mycc = c.red;
    end;
end;

function plotclrpatch(c)
figure(1);clf;
x0 = 13;y0 = 0;col = 0;
for i = 1:length(c.C),
    if rem(i,5)==1,
        col= col+1;
    end;
    row = rem(i,5);if row == 0, row = 5; end;
    xo = (col-1)*3;
    yo = (row-1)*3;
    x = xo+[x0 x0 x0+2 x0+2];
    y = yo+[y0 y0+2 y0+2 y0];
    patch(x,y,c.C{i},'EdgeColor',c.C{i});
    text(x(1)-0.2,y(1),['Colour index ' num2str(i)],'Rotation',90,'FontWeight','bold');
end;

for i = 1:length(c.grey),
    if rem(i,5)==1,
        col= col+1;
    end;
    row = rem(i,5);if row == 0, row = 5; end;
    xo = (col-1)*3;
    yo = (row-1)*3;
    x = xo+[x0 x0 x0+2 x0+2];
    y = yo+[y0 y0+2 y0+2 y0];
    patch(x,y,c.grey{i},'EdgeColor',c.grey{i});
    text(x(1)-0.2,y(1),['Grey index ' num2str(i)],'Rotation',90,'FontWeight','bold');
    if i == 1,
        text(x(1)+1,y(1)+1,'ABB','FontWeight','bold','HorizontalAlignment','center');
    end;
end;
tc{1} = c.red;
tc{2} = c.white;
tc{3} = c.black;
for i = 1:length(tc),
    if rem(i,5)==1,
        col= col+1;
    end;
    row = rem(i,5);if row == 0, row = 5; end;
    xo = (col-1)*3;
    yo = (row-1)*3;
    x = xo+[x0 x0 x0+2 x0+2];
    y = yo+[y0 y0+2 y0+2 y0];
    if i == 2, ec = c.black; else ec = tc{i};end;
    patch(x,y,tc{i},'EdgeColor',ec);
%     text(x(1)-0.2,y(1),['Grey index ' num2str(i)],'Rotation',90,'FontWeight','bold');
    if i == 1,
        text(x(1)-0.2,y(1),'''red''','Rotation',90,'FontWeight','bold');
        text(x(1)+1,y(1)+1,'ABB','FontWeight','bold','HorizontalAlignment','center');
    end;
    if i == 2,
        text(x(1)-0.2,y(1),'''white''','Rotation',90,'FontWeight','bold');
    end;
    if i == 3,
        text(x(1)-0.2,y(1),'''black''','Rotation',90,'FontWeight','bold');
    end;
end;

axis([12 43 -1 15]);
set(gca,'YTick',[]);
set(gca,'XTick',[]);
box on;


function clr = defclr
% abb blue colors
clr.C{1} = [0 75 122]/255;  % ABB_B_Head-W
clr.C{2} = [0 118 183]/255;  % ABB_B_Sub-W
clr.C{3} = [0 150 234]/255;  % ABB_B_Head-B
clr.C{4} = [91 216 255]/255;  % ABB_B_Sub-B
clr.C{5} = [172 236 255]/255;  % ABB_B_lightest

%Blue - Green colors
clr.C{6} = [0 85 95]/255;  % ABB_BG_Head-W
clr.C{7} = [0 121 135]/255;  % ABB_BG_Sub-W
clr.C{8} = [0 172 182]/255;  % ABB_BG_Head-B
clr.C{9} = [57 232 218]/255;  % ABB_BG_Sub-B
clr.C{10} = [151 255 239]/255;  % ABB_BG_lightest

%Green colors
clr.C{11} = [0 85 20]/255;  % ABB_G_Head-W
clr.C{12} = [2 130 8]/255;  % ABB_G_Sub-W
clr.C{13} = [58 178 0]/255;  % ABB_G_Head-B
clr.C{14} = [152 219 56]/255;  % ABB_G_Sub-B
clr.C{15} = [199 252 166]/255;  % ABB_G_lightest

%Green - Yellow colors
clr.C{16} = [78 86 29]/255;  % ABB_GY_Head-W
clr.C{17} = [106 116 0]/255;  % ABB_GY_Sub-W
clr.C{18} = [176 176 35]/255;  % ABB_GY_Head-B
clr.C{19} = [198 213 94]/255;  % ABB_GY_Sub-B
clr.C{20} = [224 235 146]/255;  % ABB_GY_lightest

% Yellow colors
clr.C{21} = [175 120 1]/255;  % ABB_Y_Head-W
clr.C{22} = [221 160 10]/255;  % ABB_Y_Sub-W
clr.C{23} = [238 206 35]/255;  % ABB_Y_Head-B
clr.C{24} = [244 238 127]/255;  % ABB_Y_Sub-B
clr.C{25} = [250 243 174]/255;  % ABB_Y_lightest

% Orange colors
clr.C{26} = [154 40 1]/255;  % ABB_O_Head-W
clr.C{27} = [191 69 0]/255;  % ABB_O_Sub-W
clr.C{28} = [255 108 0]/255;  % ABB_O_Head-B
clr.C{29} = [253 172 37]/255;  % ABB_O_Sub-B
clr.C{30} = [255 219 133]/255;  % ABB_O_lightest

% Orange - Violet colors
clr.C{31} = [109 21 65]/255;  % ABB_OV_Head-W
clr.C{32} = [176 74 116]/255;  % ABB_OV_Sub-W
clr.C{33} = [218 61 99]/255;  % ABB_OV_Head-B
clr.C{34} = [230 155 180]/255;  % ABB_OV_Sub-B
clr.C{35} = [241 190 197]/255;  % ABB_OV_lightest

% Violet colors
clr.C{36} = [73 56 129]/255;  % ABB_V_Head-W
clr.C{37} = [110 90 161]/255;  % ABB_V_Sub-W
clr.C{38} = [152 104 239]/255;  % ABB_V_Head-B
clr.C{39} = [180 160 232]/255;  % ABB_V_Sub-B
clr.C{40} = [216 204 248]/255;  % ABB_V_lightest

clr.white = [1 1 1];
clr.black = [0 0 0];
clr.grey{1} = [134 134 134]/255; % ABB grey
clr.grey{2}= [0.2 0.2 0.2];
clr.grey{3}= [0.4 0.4 0.4];
clr.grey{4}= [0.7 0.7 0.7];
clr.grey{5}= [0.9 0.9 0.9];
clr.red = [255 0 15]/255; % ABB red
