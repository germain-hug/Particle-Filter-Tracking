function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;

%load in the required data
load('PracticalData','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 
plot(pts1b(1,:),pts1b(2,:),'m.');
figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 
figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 
%calculate homography from pts1 to pts2, and from pts1 to pts3
H_12 = calcBestHomography(pts1, pts2);
H_13 = calcBestHomography(pts1b, pts3);

%****TO DO**** 
%for every pixel in image 1
    %transform this pixel position with your homography to find where it 
    %is in the coordinates of image 2
    %if it the transformed position is within the boundary of image 2 then 
        %copy pixel colour from image 2 pixel to current position in image 1 
        %draw new image1 (use drawnow to force it to draw)
    %end
%end;

[w1,h1,~] = size(im1);
[w2,h2,~] = size(im2);
[w3,h3,~] = size(im3);

% Compute vectors of position indices
[X,~] = meshgrid(1:h1,1:w1);
idx_i = repmat(1:w1,[1,h1]);
idx_j = reshape(X,[1,w1*h1]);

% Transform pixel position with homography 
tr_pos = H_12*[idx_j;idx_i;ones(1,size(im1,1)*size(im1,2))];

% Normalize to get homogeneous coordinates
tr_pos(1:2,:) = floor(tr_pos(1:2,:)./repmat(tr_pos(3,:),[2,1]));

for j=1:h1
    for i=1:w1
        new_x = tr_pos(2,i+(j-1)*w1);
        new_y = tr_pos(1,i+(j-1)*w1);
        % Check if transformed position is within im2 borders
        if(new_x > 0 && new_x < w2 ...
        && new_y > 0 && new_y < h2)
            im1(i,j,:) = im2(new_x,new_y,:);
        end
    end
end


%show images and points

%****TO DO****
%repeat the above process mapping image 3 to image 1.

% Transform pixel position with homography 
tr_pos = H_13*[idx_j;idx_i;ones(1,size(im1,1)*size(im1,2))];

% Normalize to get homogeneous coordinates
tr_pos(1:2,:) = floor(tr_pos(1:2,:)./repmat(tr_pos(3,:),[2,1]));

for j=1:h1
    for i=1:w1
        new_x = tr_pos(2,i+(j-1)*w1);
        new_y = tr_pos(1,i+(j-1)*w1);
        % Check if transformed position is within im2 borders
        if(new_x > 0 && new_x < w3 ...
        && new_y > 0 && new_y < h3)
            im1(i,j,:) = im3(new_x,new_y,:);
        end
    end
end




figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;



%==========================================================================
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

