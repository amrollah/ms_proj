function displayClouds(mask, bboxes)
 mask = objectAnnotation(mask,'Rectangle',bboxes,'g',3);
 frame = cv.drawContours(img,contours,'Color',[255 255 0],'Thickness',3);
    % convert the frame and the mask to uint8 RGB
    % frame = im2uint8(frame);
    % mask = uint8(repmat(mask, [1, 1, 3])) .* 255;
    % display the objects. If an object has not been detected
    % in this frame, display its predicted bounding box.
    if ~isempty(reliableTracks)
        % get bounding boxes
        bboxes = cat(1, reliableTracks.bbox);
        realKals = cat(1,reliableTracks.kalmanFilterProj);
        %!!!!! pay attention !!!!!
        realProjs = cat(1,realKals.x);
        ids = cat(1,reliableTracks.id);
        % get ids
        %                 ids = int32([reliableTracks(:).id]);
        %                 xspeeds = [reliableTracks(:).xspeed];
        %                 yspeeds = [reliableTracks(:).yspeed];

        % create labels for objects indicating the ones for
        % which we display the predicted rather than the actual
        % location
        %                 labels = cellstr(int2str(ids'));
        %                 xlspeeds = cellstr(num2str(xspeeds',2));
        %                 ylspeeds = cellstr(num2str(yspeeds',2));

        %predictedTrackInds = ...
        %    [reliableTracks(:).consecutiveInvisibleCount] > 0;
        %isPredicted = cell(size(labels));
        %isPredicted(predictedTrackInds) = {' predicted'};
        %space = {', '};
        %                 labels = strcat(labels, ', x: ', xlspeeds, ', y: ', ylspeeds);

        predictedTrackInds = ...
            [reliableTracks(:).consecutiveInvisibleCount] > 0;

        newTrackInds = ...
            [reliableTracks(:).consecutiveInvisibleCount] <= 0 & [reliableTracks(:).age] <= 5;
        trackedTrackInds = ...
            [reliableTracks(:).consecutiveInvisibleCount] <= 0 & [reliableTracks(:).age] > 5;
        if(conditionShowResults)
            if (sum(predictedTrackInds)>0)
                bboxColor = 'r';
                % draw on the mask
                mask = objectAnnotation(mask, 'Rectangle', bboxes(predictedTrackInds,:),bboxColor,3);
                %mask = objectAnnotation(mask, 'Text', [ids(predictedTrackInds,:),bboxes(predictedTrackInds,1:2)],bboxColor,3);
            end
            if (sum(trackedTrackInds)>0)
                bboxColor = 'y';
                % draw on the mask
                mask = objectAnnotation(mask, 'Rectangle', bboxes(trackedTrackInds,:),bboxColor,3);
                %mask = objectAnnotation(mask, 'Text', [ids(trackedTrackInds,:),bboxes(trackedTrackInds,1:2)],bboxColor,3);
            end
            if (sum(newTrackInds)>0)
                bboxColor = 'g';
                % draw on the mask
                mask = objectAnnotation(mask, 'Rectangle', bboxes(newTrackInds,:),bboxColor,3);
                %mask = objectAnnotation(mask, 'Text', [ids(newTrackInds,:),bboxes(newTrackInds,1:2)],bboxColor,3);
            end
        end
        % draw on the mask
        %mask = objectAnnotation(mask, 'Rectangle', bboxes,bboxColor,4);

        % draw on the frame
        frame = cv.drawContours(frame,contours,'Color',[255 255 0],'Thickness',3);
        predFrame = cv.drawContours(predFrame,contours,'Color',[0 255 0],'Thickness',3);
    end
    allPredContours = [predCloudContours,predShadowContours];
    colors = [repmat('b',size(predCloudContours)),repmat('g',size(predShadowContours))];
    % draw on the mask

    %mask = objectAnnotation(mask, 'Rectangle', bboxes);

    % draw on the frame
    %frame = cv.drawContours(frame,contours,'Color',[255 255 0]);
end