function EA_messung(file,N,EA,phi,params,fCallback,fPrepare)

global messung_i messung_file messung_xopt messung_FCount
global messung_M messung_t messung_T messung_FitStat

if nargin<6, fCallback=[]; end;
if nargin<7, fPrepare=[]; end;

if isnumeric(file), file = ['m' num2str(file)];end;
messung_file = file;
messung_xopt = {};
messung_M = [];
messung_T = [];
messung_t = [];
messung_FitStat = [];
messung_FCount = [];

if isempty(EA)
  EA = 'GA';
end;
switch lower(EA)
case 'ga'
  fEA = @GaCore;
  if isempty(fCallback) fCallback = @GA_messung_callback; end;
case 'es'
  fEA = @EsCore;
  if isempty(fCallback) fCallback = @ES_messung_callback; end;
case 'cmaes'
  fEA = @CmaesCore;
  if isempty(fCallback) fCallback = @CMAES_messung_callback; end;
otherwise error('Unknown EA');
end;

interface.DisplayGen = 0;
interface.DisplayFct = fCallback;

rand('state',sum(100*clock));
randn('state',sum(100*clock));

for messung_i=1:N  
  disp(' ');disp('******************');
  disp(['Messung Nr. ' num2str(messung_i)]);
  disp('******************');disp(' ');
  if ~isempty(fPrepare), feval(fPrepare); end;
  tic;feval(fEA, phi, params, interface);
end;
