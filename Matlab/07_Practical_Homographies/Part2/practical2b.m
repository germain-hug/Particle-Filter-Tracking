function practical2b
close all
clear all
%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
    308.9825  236.7646  255.4416  340.7335   281.5895];


%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
    50 -50 -50  50 0;...
    0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
    0    640  240;
    0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',10);

%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.
TEst = estimatePlanePose(xImCart,XCart,K);



%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
    50 -50 -50  50  50 -50 -50  50;...
    0   0   0   0 -100 -100 -100 -100];

%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points

% ---> projected points = K_Hom * TEst * XWireFrameCart_Hom
proj = [K, zeros(3,1)]*TEst*[XWireFrameCart;ones(1,size(XWireFrameCart,2))];
proj = proj(1:3,:)./repmat(proj(3,:),3,1);

% % Plotting:


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

%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?
% 5 points = over-constrained !


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

