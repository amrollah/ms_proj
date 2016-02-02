function CMAES_messung_callback

global ES_FitStat ES_Optimum ES_xOpt ES_tAct
global messung_i messung_file messung_xopt
global messung_M messung_t messung_T messung_FitStat

messung_T(messung_i) = toc;
messung_t(messung_i) = ES_tAct;
messung_M(messung_i) = ES_Optimum;
messung_FitStat(:,1:ES_tAct,messung_i) = ES_FitStat;
messung_xopt{messung_i} = ES_xOpt;

save(messung_file,'messung_M','messung_T','messung_i',...
  'messung_t','messung_FitStat','messung_xopt');
%save([messung_file '_status'],'messung_i','-ascii');
