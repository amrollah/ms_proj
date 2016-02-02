function showdbg_ppx(im)
[yy,xx]=ndgrid(im.prg_yy,im.prg_xx);
h = plot(im.ppx(:,2),im.ppx(:,1),'b.');
daspect([1 1 1]);
set(h,'ButtonDownFcn',@showdbg_ppx_click)

  function showdbg_ppx_click(varargin)
    xlim = get(gca,'xlim'); ylim = get(gca,'ylim');
    p = get(gca,'currentpoint');
    [~,j]=min(sum(bsxfun(@minus,p(1,[2 1]),im.ppx).^2,2));
    %disp([num2str(p(1,[2 1])) ' ' num2str(j)]);
    b = find(im.prg_j==j);
    h = plot(im.ppx(:,2),im.ppx(:,1),'b.',xx(b),yy(b),'ro');
    set(h,'ButtonDownFcn',@showdbg_ppx_click)
    hold on;
    if length(b)>100
      b = b(randperm(length(b)));
      b = b(1:100);
    end
    for i=1:length(b)
      [~,j] = min(sum(bsxfun(@minus,[yy(b(i)) xx(b(i))],im.ppx).^2,2));
      h2 = plot([xx(b(i)) im.ppx(j,2)],[yy(b(i)) im.ppx(j,1)],'k:');
      set(h2,'ButtonDownFcn',@showdbg_ppx_click)
    end
    h2 = plot(im.ppx(j,2),im.ppx(j,1),'gx','linewidth',2,'markersize',12);
    set(h2,'ButtonDownFcn',@showdbg_ppx_click)
    hold off;
    daspect([1 1 1]);
    set(gca,'xlim',xlim); set(gca,'ylim',ylim);
  end

end

