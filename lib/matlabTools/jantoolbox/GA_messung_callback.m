function GA_messung_callback

global GA_FitStat GA_Optimum GA_xOpt GA_tAct GA_FCount
global messung_i messung_file messung_xopt messung_FCount
global messung_M messung_t messung_T messung_FitStat

messung_T(messung_i) = toc;
messung_t(messung_i) = GA_tAct;
messung_M(:,:,messung_i) = GA_Optimum;
messung_FitStat(:,1:GA_tAct,:,:,messung_i) = GA_FitStat;
messung_xopt{messung_i} = GA_xOpt;
messung_FCount(messung_i,1:GA_tAct) = GA_FCount;

save(messung_file,'messung_M','messung_T','messung_i',...
  'messung_t','messung_FitStat','messung_xopt','messung_FCount');
%save([messung_file '_status'],'messung_i','-ascii');
