function messreihen_durchfuehren

% set random seeds first
rand('state',sum(100*clock));
randn('state',sum(100*clock));

%pause(30*rand)
prozessnr = round(rem(sum(clock),10)*100000);
nr = 0;

% Inf does not work with fscanf and fprintf, which have to be usedwith mcc !
MyInf = 9999999;

while nr<MyInf
  if exist('messreihe_stop')
    stop = load('messreihe_stop');
    if ~isempty(stop) & any(stop(:)==prozessnr), break;end;
  end;
  if exist('messreihe_status')
    fid = fopen('messreihe_status','rt');
    status = fscanf(fid,'%i %i',[2 inf])'; % status = load('messreihe_status');
    fclose(fid);
    if (length(status)>1)
      nr = status(end,1)+1;
      status = [status; nr prozessnr];
    else
      status = [1 prozessnr];
      nr = 1;
    end;
  else
    status = [1 prozessnr];
    nr = 1;
  end;
  while 1
    if nr>=MyInf, break; end;
    fid = fopen('messreihe_status','wt');
    fprintf(fid,'%i %i\n',status'); % save('messreihe_status','status','-ascii')
    fclose(fid);
    pause(2*rand+1);
    fid = fopen('messreihe_status','rt');
    status = fscanf(fid,'%i %i',[2 inf])'; % status = load('messreihe_status');
    fclose(fid);
    if status(find(status(:,1)==nr),2)==prozessnr, break; end;
    nr = status(end,1)+1;
    status = [status; nr prozessnr];
  end;

  if nr<MyInf
    disp(' ');
    disp('********************');
    disp(['Messreihe Nr. ' num2str(nr)]);
    disp('********************');
    disp(' ');
    if mess(nr)==0
      fid = fopen('messreihe_status','rt');
      status = fscanf(fid,'%i %i',[2 inf])';
      fclose(fid);
      status = [status; MyInf prozessnr];
      fid = fopen('messreihe_status','wt');
      fprintf(fid,'%i %i\n',status'); % save('messreihe_status','status','-ascii')
      fclose(fid);
    end;
  end;
end;
