function generateSynopsisTable(directory, outputFile, rootDir)
%% Create an html table of files and one line descriptions for a directory
% PMTKneedsMatlab 
%%

% This file is from matlabtools.googlecode.com

if nargin < 2, outputFile = ''; end
files        = filelist(directory, '*.m', true);
files        = files(sortidx(lower(files))); 
if isempty(files), return; end
descriptions = colvec(cellfuncell(@helpline, files)); 
fnames       = colvec(cellfuncell(@(c)c(length(directory)+2:end-2), files)); 
flinks       = cellfuncell(@(c)c(length(rootDir)+2:end), files); 
ftext        = cellfun(@googleCodeLink, flinks, fnames, 'uniformoutput', false);
tagstrings   = cellfuncell(@(f)htmlTagString(tagfinder(f)), files); 
color        = '#990000';
header       = formatHtmlText({
    '<font align="left" style="color:%s"><h2>Listing: %s</h2></font>'
    ''
    'Revision Date: %s'
    ''
    'Auto-generated by %s'
    ''
    ''
    ''
    'Some files are tagged as follows:'
    ''
    '%s'
    '<br><br><br><br><br><br><br><br><br><br>'
    ''}, color, directory(length(rootDir())+1:end), date, mfilename, htmlTagKey('left')); 
   
htmlTable('data', [ftext, descriptions, tagstrings], ...
    'colNames'          , {'FILE', 'DESCRIPTION', 'TAGS'}, ...
    'header'            , header, ...
    'colNameColors'     , {color, color, color},...
    'dataAlign'         , 'left' , ...
    'doshow'            , nargin < 2, ...
    'dosave'            , nargin > 1, ...
    'filename'          , outputFile);
end
