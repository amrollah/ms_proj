function messreihen_durchfuehren

prozessnr = round(rem(sum(clock),10)*100000);
nr = 0;

rand('state',sum(100*clock));
randn('state',sum(100*clock));

while nr<Inf
  if exist('messreihe_stop')
    stop = load('messreihe_stop');
    if ~isempty(stop) & any(stop(:)==prozessnr), break;end;
  end;
  if exist('messreihe_status')
    status = load('messreihe_status');
    nr = status(end,1)+1;
    status = [status; nr prozessnr];
  else
    status = [1 prozessnr];
    nr = 1;
  end;
  while 1
    if nr==Inf, break; end;
    save('messreihe_status','status','-ascii')
    pause(2);
    status = load('messreihe_status');
    if status(find(status(:,1)==nr),2)==prozessnr, break; end;
    nr = status(end,1)+1;
    status = [status; nr prozessnr];
  end;
  
  if nr<Inf
    disp(' ');
    disp('********************');
    disp(['Messreihe Nr. ' num2str(nr)]);
    disp('********************');
    disp(' ');
    if mess(nr)==0
      status = [load('messreihe_status'); Inf prozessnr];
      save('messreihe_status','status','-ascii');
    end;
  end;
end;
