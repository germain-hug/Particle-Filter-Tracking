function r=HW2_Practical9c( templateMetaFileName )
%
% You'll end up calling this function four times (once each with 'll',
% 'lr', 'ul', 'ur'). For each time you'll get an estimated track (2D
% locations over time) of corner of the pattern moved. Note that we're
% saving off the single MAP estimate at each frame, but the particles in
% the Particle Filter are attempting to keep track of the whole posterior
% distribution.
%
% Most TO DO parts below should be copy-paste from Practical9a and
% Practical9b.
%
% Some implementation trivia (i.e. you shouldn't have to worry about
% these changes):
% MeasurePatchSimilarityHere.m and PatchDistPicker.m have been modifid from
% their previous versions to cope with these particular gray-scale images.
% If you want to run ExhaustiveTemplateSearchDemo.m, you'll want to run
% with the older versions of those two files (available on Moodle).

LoadVideoFrames

% Load template and starting position ('pos').
%load ll;
%load lr;
%load ul;
%load ur;
load( templateMetaFileName );


numParticles = 100;
weight_of_samples = ones(numParticles,1);

% TO DO: normalize the weights (may be trivial this time)
weight_of_samples = weight_of_samples/sum(weight_of_samples);

% Initialize which samples from "last time" we want to propagate: all of
% them!:
samples_to_propagate = [1:numParticles]';


% ============================
% NOT A TO DO: You don't need to change the code below, but eventually you may
% want to vary the number of Dims.
numDims_w = 2;

% NOT A TO DO: You don't need to change anything here but please make sure you
% MATLAB console to find more info for the buildin function.
%
% Here we randomly initialize some particles throughout the space of w.
% The positions of such particles are quite close to the know initial
% position:
particles_old = repmat([minY minX], numParticles, 1 ) + 10*rand( numParticles, numDims_w );
% ============================


hImg = figure;
hSamps = figure;

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Was our "Loop from here" point in Part a.

% Initialize a temporary array r to store the per-frame MAP estimate of
% w. This is what we'll return in the end.
r = zeros(numFrames, numDims_w);

for( iTime = 1:numFrames )
    
    % TO DO: compute the cumulative sume of the weights. You could refer to
    % the MATLAB buitin function 'cumsum'.
    cum_hist_of_weights = cumsum(weight_of_samples);
    
    
    % ==============================================================
    % Resample the old distribution at time t-1 and select samples favoring
    % those with a higher posterior probability.
    % ==============================================================
    samples_to_propagate = zeros(numParticles,1);
    
    % Pick random thresholds (numParticles X 1 matrix) in the cumulative
    % probability's range [0,1]:
    some_threshes = rand(numParticles,1);
    
    
    % For each random threshold, find which sample in the ordered set is
    % the first one to push the cumulative probability above that
    % threshold, e.g. if the cumulative histogram goes from 0.23 to 0.26
    % between the 17th and 18th samples in the old distribution, and the
    % threshold is 0.234, then we'll want to propagate the 18th sample's w
    % (i.e. particle #18).
    
    for sampNum = 1:numParticles
        thresh = some_threshes(sampNum);
        for index = 1:numParticles
            if( cum_hist_of_weights(index) > thresh )
                break;
            end;
        end;
        samples_to_propagate(sampNum) = index;
    end;
    % Note: it is fine even if some of the old particles get picked repeatedly, while
    % others don't get picked at all.
    
    
    % =================================================
    % Visualize
    % =================================================
    set(0,'CurrentFigure',hSamps)
    set(gcf,'Position',[600 125 640 480]);
    set(gcf,'Color',[1 1 1]);
    set(gca,'Box','Off');
    title('Cumulative histogram of probabilities for sorted list of particles');
    plot(zeros(numParticles,1), some_threshes,'b.','MarkerSize',15);
    hold on;
    plot([1:numParticles], cum_hist_of_weights, 'ro-', 'MarkerSize',3);
    which_sample_ids = unique(samples_to_propagate);
    how_many_of_each = histc(samples_to_propagate, unique(samples_to_propagate));
    for( k=1:size(which_sample_ids,1) )
        plot( which_sample_ids(k), 0, 'bo-', 'MarkerSize', 3 * how_many_of_each(k) )
    end;
    xlabel(sprintf( 'Indeces of all available samples, with larger blue circles for frequently re-sampled particles\n(Iteration %d)', iTime));
    ylabel('Cumulative probability');
    hold off
    xlim([0 numParticles]);ylim([0 1]);
    drawnow
    % =================================================
    % =================================================
    
    
    % Predict where the particles we plucked from the old distribution of
    % state-space will go in the next time-step. This means we have to apply
    % the motion model to each old sample.
    
    particles_new = zeros( size(particles_old) );
    for particleNum = 1:numParticles
        
        % TO DO: Incorporate some noise, e.g. Gaussian noise with std 10,
        % into the current location (particles_old), to give a Brownian
        % motion model. -> we randomize 2D position
        noise = randn(1,size(particles_old,2))*10;
        particles_new(particleNum,:) = particles_old( samples_to_propagate(particleNum),: ) + noise;
        particles_new(particleNum,:) = round(  particles_new(particleNum,:)  ); % Round the particles_new to simplify Likelihood evaluation.
    end;
    
    
    % TO DO (optional): change the motion model above to have different
    % degrees of freedom, and optionally completely different motion models.
    % See Extra Credit for more instructions.
    
    Im2   = double( Imgs{iTime} );
    Edge2 = edge(Im2);

    % From here we incorporate the measurement for the new state (time t):
    % The new particles at predicted locations in state-space
    % for time t, are missing their weights: how well does each particle
    % explain its observation x_t?
    for particleNum = 1:numParticles
        
        % Convert the particle from state-space w to measurement-space x:
        % Note: Trivial here because both are in 2D space of image
        % coordinates
        
        % Within the loop, we evaluate the likelihood of each particle:
        particle = particles_new(particleNum,:);
        % Check that the predicted location is a place we can really evaluate
        % the likelihood.
        s = size(pixelsTemplate);
        inFrame = particle(1) >= 1.0   &&  particle(1)+ s(1) <= imgHeight && ...
            particle(2) >= 1.0   &&  particle(2) + s(2) <= imgWidth;
        if( inFrame )
            minX = particles_new(particleNum,2);
            minY = particles_new(particleNum,1);
            
            % A stand-alone demo (ExhaustiveTemplateSearchDemo.m) is
            % provided to illustrate the functionality of our likelihood
            % function:
            % MeasurePatchSimilarityHere(im,template,x,y)
            % im : image to search
            % template : the given template
            % x, y : the upper left corner in the image to place the
            % template for similarity evaluation
            weight_of_samples(particleNum) = ...
                MeasurePatchSimilarityHere( Im2, pixelsTemplate, minY, minX );
            if(minX > 0   &&  minX <= size(Edge2,1) && minY > 0   &&  minY <= size(Edge2,2) )
                weight_of_samples(particleNum) = weight_of_samples(particleNum)...
                    + 0.5*Edge2(minX,minY);
            end
        else
            weight_of_samples(particleNum) = 0.0;
        end;
        
    end;
    
    % TO DO: normalize the weights
    weight_of_samples = weight_of_samples/sum(weight_of_samples);
    
    % TO DO: Compute the coordinate of the "best" (i.e. MAP) location by
    % computing the weighted average of all the particles:
    weightedAve = sum(repmat(weight_of_samples,1,2).*particles_new);
    
    middleOfTrackedPatch = weightedAve + patchOffset;
    r(iTime,:) = middleOfTrackedPatch; % Return the MAP of middle position
    % of a patch. That 'middle' is the patch's interest point in the
    % current frame. The coordinate we've been using everywhere is the
    % upper-left corner of each patch. We're actually most interested
    % in the corner of the square. In this case, a patchOffset that goes
    % along with each pixelsTemplate, is given to reset the interest point
    % back to the middle of the patch.
    
    set(0,'CurrentFigure',hImg)
    imagesc(Im2/255)
    colormap(gray);
    set(gcf,'Position',[23 125 640 480]);
    set(gcf,'Color',[1 1 1]);
    title(sprintf( 'Particles projected to measurement-space\n(Time %d)', iTime));
    hold on
    plot(particles_new(:,2), particles_new(:,1), 'rx')
    plot(middleOfTrackedPatch(2), middleOfTrackedPatch(1), 'bo')
    hold off
    drawnow
    
    % Optional code to save out figure:
    %     pngFileName = sprintf( '%s_%.5d.png', 'myOutput', iTime );
    %     saveas(gcf, pngFileName, 'png');
    %     print( gcf, '-dpng', '-r80', pngFileName ); % Gives 640x480 (small) figure
    
    
    % Now we're done updating the state for time t.
    % For Condensation, just clean up and prepare for the next round of
    % predictions and measurements:
    particles_old = particles_new;
    clear particles_new;
    
end; % End of for loop over each frame in the sequence
