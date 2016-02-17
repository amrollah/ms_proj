function check_segmentation(interactive)
if nargin<1, interactive=false; end
if interactive
  show1 = @show_interactive;
else
  show1 = @show_simple;
end
figcount = 0;

s = vmlSeq('cavriglia','2015_07_23');
if interactive, figure(111);s.plotpow; end
s.loadframe(5318); feval(show1);
s.loadframe(5418); feval(show1);
s.loadframe(5500); feval(show1);

s = vmlSeq('cavriglia','2015_07_29');
if interactive, figure(111);s.plotpow; end
s.loadframe(4001); feval(show1);

s = vmlSeq('cavriglia','2015_07_31');
if interactive, figure(111);s.plotpow; end
s.loadframe(4820); feval(show1);
s.loadframe(3820); feval(show1);
s.loadframe(2820); feval(show1);
s.loadframe(1820); feval(show1);
s.loadframe(820); feval(show1);
s.loadframe(1000); feval(show1);
s.loadframe(2000); feval(show1);
s.loadframe(3000); feval(show1);
s.loadframe(4000); feval(show1);
s.loadframe(5000); feval(show1);
s.loadframe(6000); feval(show1);

s = vmlSeq('cavriglia','2015_08_04');
if interactive, figure(111);s.plotpow; end
s.loadframe(1000); feval(show1);
s.loadframe(2000); feval(show1);
s.loadframe(3000); feval(show1);
s.loadframe(4000); feval(show1);
s.loadframe(5000); feval(show1);
s.loadframe(6000); feval(show1);

s = vmlSeq('cavriglia','2015_08_09');
if interactive, figure(111);s.plotpow; end
s.loadframe(4401); feval(show1);
s.loadframe(4402); feval(show1);
s.loadframe(4539); feval(show1);
s.loadframe(4570); feval(show1);

s = vmlSeq('cavriglia','2015_08_16');
if interactive, figure(111);s.plotpow; end
s.loadframe(5000); feval(show1);
s.loadframe(5001); feval(show1);
s.loadframe(5002); feval(show1);
s.loadframe(5222); feval(show1);
s.loadframe(6500); feval(show1);
s.loadframe(6501); feval(show1); %presently does not work well!

% s = vmlSeq('cavriglia','2015_08_20');
% if interactive, figure(111);s.plotpow; end
% s.loadframe(5538); feval(show1);

s = vmlSeq('cavriglia','2016_01_11');
if interactive, figure(111);s.plotpow; end
s.loadframe(1); feval(show1);
s.loadframe(300); feval(show1);
s.loadframe(1607); feval(show1);

  function show_interactive
    s.ShowSeg; pause;
  end

  function show_simple
    figcount = figcount+1;
    figure(figcount); clf;
    s.showseg;
  end
end
