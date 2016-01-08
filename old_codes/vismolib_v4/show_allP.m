for i=1:length(allP)
  if ~isempty(allP{i})
    figure(100);clf;
    plot(allP{i}(:,1),allP{i}(:,2));
    datetickzoom;
    title([num2str(i) ' ' folders{i} ' ' num2str(Pncol(i))],'interpreter','none');
    pause;
  end
end
