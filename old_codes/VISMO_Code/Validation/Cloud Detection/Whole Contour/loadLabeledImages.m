HOMEANNOTATIONS = '..\..\Utils\LabeledImages\Annotations'; % set this to the folder in which we have the labels\annotations
HOMEIMAGES = '..\..\Utils\LabeledImages\Images'; % set this to the folder in which we have the images
folderlist = {'users\burakzeydan\cloud_images'};
%LMinstall(folderlist,HOMEIMAGES,HOMEANNOTATIONS)
D = LMdatabase(HOMEANNOTATIONS,folderlist);
LMdbshowscenes(D,HOMEIMAGES);
LMplot(D,1,HOMEIMAGES);
%cloudDetection_Validation(D,1,HOMEIMAGES)