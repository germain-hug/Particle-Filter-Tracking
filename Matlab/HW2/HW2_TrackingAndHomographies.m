function HW2_TrackingAndHomographies
close all;
clear all;

LLs = HW2_Practical9c( 'll' );
LRs = HW2_Practical9c( 'lr' );
ULs = HW2_Practical9c( 'ul' );
URs = HW2_Practical9c( 'ur' );

close all;

% Load frames from the whole video into Imgs{}.
% This is really wasteful of memory, but makes subsequent rendering faster.
LoadVideoFrames

% Coordinates of the known target object (a dark square on a plane) in 3D:
XCart = [-50 -50  50  50;...
    50 -50 -50  50;...
    0   0   0   0];

% These are some approximate intrinsics for this footage.
K = [640  0    320;...
    0    512  256;
    0    0    1];

% Define 3D points of wireframe object.
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
    50 -50 -50  50  50 -50 -50  50;...
    0   0   0   0 -100 -100 -100 -100];

hImg = figure;

% ================================================
for iFrame = 1:numFrames
    xImCart = [LLs(iFrame,:)' ULs(iFrame,:)' URs(iFrame,:)' LRs(iFrame,:)'];
    xImCart = circshift( xImCart, 1);
    
    % To get a frame from footage
    im = Imgs{iFrame};
    
    % Draw image and 2d points
    set(0,'CurrentFigure',hImg);
    set(gcf,'Color',[1 1 1]);
    imshow(im); axis off; axis image; hold on;
    plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',15);
    
    
    %TO DO Use your routine to calculate TEst the extrinsic matrix relating the
    %plane position to the camera position.
    T = estimatePlanePose(xImCart, XCart, K);
    
    
    
    %TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
    %through the projective camera, and drawing lines betweeen the
    %resulting 2d image points
    hold on;
    
    % TO DO: Draw a wire frame cube using data XWireFrameCart. You need to
    % 1) project the vertices of a 3D cube through the projective camera;
    % 2) draw lines betweeen the resulting 2d image points.
    % Note: CONDUCT YOUR CODE FOR DRAWING XWireFrameCart HERE
    
    % Projection ---> projected points = K_Hom * TEst * XWireFrameCart_Hom
    proj = [K, zeros(3,1)]*T*[XWireFrameCart;ones(1,size(XWireFrameCart,2))];
    proj = proj(1:3,:)./repmat(proj(3,:),3,1);
    
    % Plotting
    plot(proj(1,:),proj(2,:),'rx')
    plot([proj(1,1),proj(1,2)],[proj(2,1),proj(2,2)],'r-')
    plot([proj(1,2),proj(1,3)],[proj(2,2),proj(2,3)],'r-')
    plot([proj(1,1),proj(1,4)],[proj(2,1),proj(2,4)],'r-')
    plot([proj(1,3),proj(1,4)],[proj(2,3),proj(2,4)],'r-')
    plot([proj(1,1),proj(1,5)],[proj(2,1),proj(2,5)],'r-')
    plot([proj(1,2),proj(1,6)],[proj(2,2),proj(2,6)],'r-')
    plot([proj(1,3),proj(1,7)],[proj(2,3),proj(2,7)],'r-')
    plot([proj(1,4),proj(1,8)],[proj(2,4),proj(2,8)],'r-')
    plot([proj(1,5),proj(1,6)],[proj(2,5),proj(2,6)],'r-')
    plot([proj(1,6),proj(1,7)],[proj(2,6),proj(2,7)],'r-')
    plot([proj(1,5),proj(1,8)],[proj(2,5),proj(2,8)],'r-')
    plot([proj(1,7),proj(1,8)],[proj(2,7),proj(2,8)],'r-')

    
    
    hold off;
    drawnow;
    
%         Optional code to save out figure
         pngFileName = sprintf( '%s_%.5d.png', 'myOutput', iFrame );
         print( gcf, '-dpng', '-r80', pngFileName ); % Gives 640x480 (small) figure
    
    
end % End of loop over all frames.
% ================================================

% TO DO: QUESTIONS TO THINK ABOUT...

% Q: Do the results look realistic?
% If not then what factors do you think might be causing this


% TO DO: your routines for computing a homography and extracting a
% valid rotation and translation GO HERE. Tips:
%
% - you may define functions for T and H matrices respectively.
% - you may need to turn the points into homogeneous form before any other
% computation.
% - you may need to solve a linear system in Ah = 0 form. Write your own
% routines or using the MATLAB builtin function 'svd'.
% - you may apply the direct linear transform (DLT) algorithm to recover the
% best homography H.
% - you may explain what & why you did in the report.



% ============================================================
% ============================================================

function T = estimatePlanePose(xImCart,XCart,K)


%TO DO Convert Cartesian image points xImCart to homogeneous representation
%xImHom
xImHom = [xImCart;ones(1,size(xImCart,2))];

%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom
xCamHom = K^(-1)*xImHom;

%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
H = calcBestHomography(XCart,xCamHom);

%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD
Phi_est = H(:,1:2);
[U,~,V] = svd(Phi_est);
L = [1 0; 0 1; 0 0];
R = U*L*V';

%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns
R = [R, cross(R(:,1),R(:,2))];
%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
if(det(R)<0)
    R(:,3) = -R(:,3);
end

%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third column of H
lambda = sum(sum(H(:,1:2)./R(:,1:2)))/6;
t = H(:,3)/lambda;

%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.
if(t(3)<0)
    t=-t;
    R(:,1:2) = -R(:,1:2);
end

%assemble transformation into matrix form
T  = [R t;0 0 0 1];


% ============================================================
% ============================================================

function H = calcBestHomography(pts1Cart, pts2Cart)

pts1Cart = [pts1Cart; ones(1,size(pts1Cart,2))];
pts2Cart = [pts2Cart; ones(1,size(pts2Cart,2))];

A = zeros(2*size(pts1Cart,2),9);
for i = 1:size(pts1Cart,2)
    ui = pts1Cart(1,i);
    vi = pts1Cart(2,i);
    xi = pts2Cart(1,i);
    yi = pts2Cart(2,i);
    A(2*i-1,:) = [0,0,0,-ui,-vi,-1,yi*ui,yi*vi,yi];
    A(2*i,:) = [ui,vi,1,0,0,0,-xi*ui,-xi*vi,-xi];
end

h = solveAXEqualsZero(A);

H = reshape(h,[3,3])';

%==========================================================================
function x = solveAXEqualsZero(A)

%****TO DO **** Write this routine
[~,~,V] = svd(A);
x = V(:,end);



