function score = MeasurePatchSimilarityHere( Im2, pixelsTemplate, minY, minX )

SearchPadding = [0 0 0 0];
[cropCoordsFullI availPadW availPadN] = ...
    IndexesToSearchInFullImg(size(Im2), size(pixelsTemplate), [minY minX], SearchPadding);
subI = Im2(cropCoordsFullI(1):cropCoordsFullI(2), cropCoordsFullI(3):cropCoordsFullI(4));
% Note: The line above implicitly gives the search sub-image, subI, an
% extra row and an extra column if available.

% imgMaskCel = GetCelMaskPixels( maskSearchWhat, 1, minXs, maxXs, minYs, maxYs );

% Just use green channel for now (ie :,:,2).
% Note!! Actually, code below has been modified to assume there is NO COLOR CHANNEL!
[ncc movedSE scoreMaxs] = ...
    PatchDistPicker( pixelsTemplate(:,:), subI(:,:), ...
                     [availPadW availPadN], ...
                     'plain_nssd', ...   % builtinNCC, builtinNCCsearchInside, plain_ncc, maskedNCC, plain_nssd, masked_nssd
                     0 );
                 
if( numel(scoreMaxs) )
    score = max( 0, scoreMaxs(1) ); % Calling functions don't like negative scores.
else
    score = 0;
end;