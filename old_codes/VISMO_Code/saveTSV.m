function saveTSV(Data, fileName)

    fileID = fopen(fileName, 'w'); %warning, existing content of the same file will be deleted.
    title = 'times\tMeasured GHI\tClear Sky Model';
    fprintf(fileID, title);
    for i= 1:size(Data{1},1) %TODO: this can be speeded up though it takes less than 2 seconds now.
       if(isnan(Data{3}(i,:)))
           fprintf(fileID, '\n%8s\t%.0f\t%.0f',Data{1}(i,:), Data{2}(i,:),0);
       else
           fprintf(fileID, '\n%8s\t%.0f\t%.0f',Data{1}(i,:), Data{2}(i,:),Data{3}(i,:));
       end
    end
    fclose(fileID);
    fprintf('file print for %s is complete\n', fileName);
end