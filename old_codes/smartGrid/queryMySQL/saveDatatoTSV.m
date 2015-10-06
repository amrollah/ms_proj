function saveDatatoTSV(Data, ColumnNames, fileName)
    fileID = fopen(fileName, 'w'); %warning, existing content of the same file will be deleted.
    
    for i=1:length(ColumnNames) %write the column titles.
        if(i == length(ColumnNames))
            fprintf(fileID, '%s', char(ColumnNames{i}));
        else
            fprintf(fileID, '%s\t', char(ColumnNames{i}));
        end
    end
    fieldNames = fieldnames(Data{1});
    %NOTICE that I hard coded the number of columns, if you plan to add more
    %into the same tsv, change here.
    for i=1:length(Data) %iterate over different timestamps
        for col=1:length(ColumnNames) %iterate over columns
            if (col == 1) %time (fieldname 1)
                fprintf(fileID, '\n%s\t', datestr(Data{i}.(fieldNames{1}),31));
            elseif (col == 2) %fieldName 2
                fprintf(fileID, '%.4f\t', Data{i}.(fieldNames{2}));
            elseif (col == 3) %fieldName 3
                fprintf(fileID, '%.4f', Data{i}.(fieldNames{3}));
            end
        end
    end
    fclose(fileID);
    fprintf('file print for %s is complete\n', fileName);
end