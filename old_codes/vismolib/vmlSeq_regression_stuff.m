%% prediction and training

    function y = pred(obj, frame_nos, PBP, W1, Z1, P1)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      if nargin<3, PBP = []; end
      if nargin<4, W1 = []; end
      if nargin<5, Z1 = []; end
      if nargin<6, P1 = []; end
      if isempty(W1), 
        [jW, tpred, ~, usepow] = obj.getWidx(PBP);
        W1 = obj.W{jW}; 
      else
        tpred = PBP(1); usepow = PBP(3);
      end
      assert(length(W1)==11+usepow);
      if isempty(Z1)
        Z1 = obj.Z{obj.getZidx(tpred)};
        Z1 = Z1(:,frame_nos);
      end
      N = obj.conf.pred.ngrid;
      [xx,yy]=ndgrid((-(N-1)/2:(N-1)/2)*2/(N-1),(-(N-1)/2:(N-1)/2)*2/(N-1));
      WW = exp(-0.5*(((xx-W1(3))/W1(4)).^2+((yy-W1(5))/W1(6)).^2));
      WW = WW(:)/sum(WW(:));
      y = tanh(W1(7)*WW'*tanh(50*(Z1(1:N^2,:)-W1(1)))+...
          W1(8)*WW'*tanh(50*(Z1(N^2+1:end,:)-W1(2)))+W1(9))*W1(10)+W1(11);
      if usepow 
        if isempty(P1), P1 = obj.getP(obj.ti(frame_nos)); end
        y = y+W1(12)*P1; 
      end
    end
    
    function y1 = get_training_target(obj,frame_nos,PBP)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      if nargin<3, PBP=[]; end
      if isempty(PBP), [~,tpred,tblur] = obj.getWidx(PBP);
      else tpred = PBP(1); tblur = PBP(2); end
      if tblur==0
        y1 = obj.getP(obj.ti(frame_nos)+tpred/86400);
      else
        y1 = zeros(1,length(frame_nos));
        for i=1:length(frame_nos)
          y1(i) = sign(tblur)*max(sign(tblur)*obj.getP(...
            obj.ti(frame_nos(i))+(tpred-abs(tblur):tpred+abs(tblur))/86400));
        end
      end
    end
    
    function err = mse(obj, frame_nos, PBP, W1, y1, Z1, P1)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      if nargin<3, PBP = []; end
      if nargin<4, W1 = []; end
      if nargin<5, y1 = []; end
      if nargin<6, Z1 = []; end
      if nargin<7, P1 = []; end
      y = obj.pred(frame_nos, PBP, W1, Z1, P1);
      if isempty(y1), y1 = obj.get_training_target(frame_nos, PBP); end
      diffy = y-y1;
      err = mean(diffy(~isnan(diffy)).^2);
      %disp(err);
    end
    
    function err = mse_persist(obj, tpred, frame_nos)
      if nargin<2 || isempty(tpred), [~,tpred] = obj.getZidx([]); end
      if nargin<3 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      y = obj.getP(obj.ti(frame_nos));
      y1 = obj.getP(obj.ti(frame_nos)+tpred/86400);
      err = mean((y-y1).^2);
    end
    
    function err = trainerr(obj, frame_nos, PBP, W1, y1, Z1, P1)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      if nargin<3, PBP = []; end
      if nargin<4, W1 = []; end
      if nargin<5, y1 = []; end
      if nargin<6, Z1 = []; end
      if nargin<7, P1 = []; end
      y = obj.pred(frame_nos, PBP, W1, Z1, P1);
      if isempty(y1), y1 = obj.get_training_target(frame_nos, PBP); end
      err = y-y1;
    end
    
    function train(obj, frame_nos, PBP, do_replace)
      if nargin<2 || strcmpi(frame_nos,'all')
        frame_nos = 1:length(obj.ti); 
      end
      if nargin<3 || isempty(PBP), PBP = [obj.Z_tpred(1) 0 0]; end
      if nargin<4, do_replace = 0; end
      jW = obj.getWidx(PBP, 1);
      y1 = obj.get_training_target(frame_nos, PBP);
      Z1 = obj.Z{obj.getZidx(PBP(1))};
      Z1 = Z1(:,frame_nos);
      jdel = isnan(y1) | any(isnan(Z1),1);
      y1(jdel) = []; Z1(:,jdel) = []; frame_nos(jdel) = [];
      P1 = obj.getP(obj.ti(frame_nos));
      maxP = max(obj.P);
      wmin = [0 0 -1 0.1 -1 0.1 -1000 -1000 -1000 -2*maxP -2*maxP -2];
      wmax = [1 1 1 10 1 10 1000 1000 1000 2*maxP 2*maxP 2];
      w0 = [0.65 0.65 0 0.5 0 0.5 1 1 0 maxP/2 maxP/2 0];
      wstep = [0.1 0.1 0.1 0.1 0.1 0.1 1 1 1 maxP/10 maxP/10 0.1];
      jdeact = [];
      wmin(jdeact)=w0(jdeact); wmax(jdeact)=w0(jdeact);
      if ~PBP(end)
        wmin(end)=[]; wmax(end)=[]; w0(end)=[]; wstep(end)=[];
      end
      tic;
      if obj.conf.trainreg.use_bobyqa
        [w, ~, retcode] = bobyqa(@(w)obj.mse(frame_nos,PBP,w,y1,Z1,P1), ...
          w0, wmin, wmax, wstep, 1000);
      else
        [w,~,~,retcode] = lsqnonlin(@(w)obj.trainerr(frame_nos,PBP,w,y1,Z1,P1), ...
          w0, wmin, wmax);
      end
      rmse = sqrt(obj.mse(frame_nos,PBP,w,y1,Z1,P1));
      if obj.printlevel>=1
        disp(['training RMSE = ' num2str(rmse) ' / retcode = ' num2str(retcode) ' / elapsed = ' num2str(toc)]);
      end
      if isempty(jW) 
        jW=size(obj.W_PBP,2)+1;
        obj.W_PBP(:,jW) = PBP; obj.W_RMSE(jW) = inf; 
      end
      if rmse<obj.W_RMSE(jW) || do_replace
        obj.W{jW} = w;
        obj.W_RMSE(jW) = rmse;
      end
    end

% end prediction and training
