% LoadVideoFrames.m
% This optional script may help you get started with loading of photos.
%
% For simplicity, this example code just
% makes a few cell-arrays to hold the hard-coded filenames.

VideoDirName = '../Pattern02';

if( ~exist( VideoDirName, 'dir') )
    display('Please change current directory to the parent folder above %s', VideoDirName);
end

allJpgs = dir(fullfile( VideoDirName, '*.jpg'));
numFrames = size(allJpgs, 1);
frameNames = cell( numFrames, 1);

for iFrame = 1:numFrames
    frameNames{iFrame} = fullfile( VideoDirName, allJpgs(iFrame).name );
end

clear allJpgs VideoDirName


Imgs = cell(numFrames,1);
for (iFrame = 1:numFrames)
    Itemp = imread( frameNames{iFrame} );
    Imgs{iFrame} = Itemp(:,:,1);
end;
imgWidth = size(Imgs{1}, 2);
imgHeight = size(Imgs{1}, 1);

clear iFrame Itemp;