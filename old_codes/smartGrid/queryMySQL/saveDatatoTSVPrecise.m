function saveDatatoTSVPrecise(Data, ColumnNames, fileName)
    fileID = fopen(fileName, 'w'); %warning, existing content of the same file will be deleted.
    
    for i=1:length(ColumnNames) %write the column titles.
        if(i == length(ColumnNames))
            fprintf(fileID, '%s', char(ColumnNames{i}));
        else
            fprintf(fileID, '%s\t', char(ColumnNames{i}));
        end
    end
    fieldNames = fieldnames(Data);
    %NOTICE that I hard coded the number of columns, if you plan to add more
    %into the same tsv, change here.
    for i=1:length(Data.Time) %iterate over different timestamps
        for col=1:length(ColumnNames) %iterate over columns
            if (col == 1) %time (fieldname 1)
                fprintf(fileID, '\n%s', datestr(Data.(fieldNames{col})(i),31));
            elseif (col == 2) %fieldName 2
                if(length(Data.(fieldNames{col}))>= i)
                    fprintf(fileID, '\t%.4f', Data.(fieldNames{col})(i));
                else
                    fprintf(fileID, '\t');
                end
            elseif (col == 3) %fieldName 3
                if(length(Data.(fieldNames{col}))>= i)
                    fprintf(fileID, '\t%.4f', Data.(fieldNames{col})(i));
                else
                    fprintf(fileID, '\t');
                end
            end
        end
    end
    fclose(fileID);
%     fprintf('file print for %s is complete\n', fileName);
end