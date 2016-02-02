function xlimit = confHullGlobal_1(xd,xm,xt,type,r,BPdim)
%
% Prozedur 			confHullGlobal_1
% Aufgabe 			Auswertung einer Hüllenfunktion mit Konfidenzwert
%               Die Hülle selbst wird durch einen Abschneidealgorithmus 
%               aus dem Hyperkubus berechnet. Dabei werden die Hyperebenen
%               senkrecht zur Verbindungslinie Zentralpunkt Constraint-Punkt
%               durch den Constraint-Punkt verlaufend gelegt.
% 						
% In 	  	      xd:    Constraint-Punkte zur Bestimmung der Hülle
%               xm:    Zentralpunkt innerhalb des Hyperkubus
%               xt:    Zu testende Punkte
%               type:  Konfidenztyp 
%                      0: einfacher Abschneidealgorithmus
%                      1: Exp-Funktion
%                      2: Konfidenzwert mit Exp-Funktion
%                      3: Konfidenzwert mit Sinus-Funktion
%               r:     Radius für Konfidenzwertberechnung
%               BPdim: Die Komponenten 1:BPdim entsprechen der Betriebsebene
% Out           D:     Abstand von den Facetten (angepasst mit Konfidenzwert)
%               I:     Indizes der Constraint-Punkte
%
% Seiteneffekte 
% Author   			Kosmas Knödler
% Version  			1.0.0
% Stand    			16/4/02
% Historie
% Ersterstellung: 04/2002, V 1.0.0
% Änderungen:
%
% Beschreibung:
%
% See also 

if nargin<6, BPdim = 0; end;
if nargin<5, r = 0.1; end;
if nargin<4, type = 0; end;


[ndim, nt] = size(xt);
nd = size(xd,2);
rixt = ones(1,nt);
rixd = ones(1,nd);
ridim = ones(1,ndim);

dmax = sqrt(ndim/2);

if size(xm,2)>1
  BP = xm(1:BPdim);
  BP = BP(:,ones(1,size(xd,2)));
  XD = [BP; xd(BPdim+1:end,:)];
else
  XD = xd;
end;

% Zentralpunkt xm wird neuer Nullpunkt
XD = XD - xm(:,rixd);
XT = xt - xm(:,rixt);

% Abstände der Punkte vom Zentralpunkt
XDf = sqrt(sum(XD.^2));
XTf = sqrt(sum(XT.^2));

% Normierung der Punkte
XDn = XD./XDf(ridim,:);
XTn = XT./XTf(ridim,:);

D = reshape(XTn,ndim,1,nt);
D = squeeze(sum((D(:,rixd,:)-XDn(:,:,rixt)).^2,1).^(1.5));
j = any(D<eps,1);
jq = find(~j);
if ~isempty(jq)
  D(:,jq)=1./D(:,jq);
  dsum = sum(D(:,jq),1);
  D(:,jq) = D(:,jq)./dsum(rixd,:);
end;
j = find(j);
if ~isempty(j)
  [v,jj] = min(D(:,j),[],1);
  D(:,j) = 0;
  D((j-1)*nd+jj) = 1;
end;

xdf = dmax-simcmodel(XDn,XDn,dmax-XDf,0.1);  
xdf_d = (XDf-xdf)./XDf;
xdn_korr = XDn.*(1+xdf_d(ridim,:));

if ndim>1, rr = 1./((1+4*XDf).^(1/(ndim-1))); else rr=1; end;
conf1 = simconfmodel(XTn,xdn_korr,rr);
conf2 = simconfmodel(XTn,XDn,rr);

y2 = sum(XDf(rixt,:)'.*D)+dmax*(1-conf2);
y = sum(xdf(rixt,:)'.*D);
y = y+(y2-y).*(1-conf1);

xlimit = xm(:,rixt)+y(ridim,:).*XTn;