clear;
clc;
% read image
[imname,impath]=uigetfile({'*.jpg;*.png;*.jpeg'}); %read all image extensions
im=imread([impath,'/',imname]); %get absolute path of image
load db; %load database
warning('off'); %ignore unnecessary warnings in output

%-------------------------------- Image Pre-processing -----------------------------%
%resize image
im=imresize(im,[128 128]); %resize the image

%remove noise;
%seperate color into red, green, blue
red=im(:,:,1);
green=im(:,:,2);
blue=im(:,:,3);

%remove noise from each color system
red=medfilt2(red);
green=medfilt2(green);
blue=medfilt2(blue);

%restore channels
channels(:,:,1)=red;
channels(:,:,2)=green;
channels(:,:,3)=blue;

figure
subplot(121)
imshow(im);
title('Resized original image');
subplot(122)
imshow(channels);
title('Resized filtered image');

%feature extraction
attrb=featureExtraction(channels);

l=length(currency);

%uses euclidean weight distance to find the highest matching currency
for c=1:l
    distance(c)=dist(attrb', currency(c).feature);  
end

[value,index]= min(distance);

if value
    nameOfCurrency = currency(index).name;
    fprintf('Recognized currency is : ');
    disp(nameOfCurrency)
else
    disp('Currency not recognized');
end
 


function attrb=featureExtraction(value)
 %color feature
 color=colorDetection(value);
 
 %edge feature
 edge=edgeDetection(value);
 
 %texture feature
 %glcm-gray level co-occurrence matrix
 glcm=graycomatrix(rgb2gray(value));
 gray_level=glcm(:);
 
 attrb=[color; edge; gray_level];
 
end




function extract = colorDetection(img)
%Obtain color moments for mean, variance and skewness in each channel from the RGB image

%convert rgb to xyz colorspace
xyzT = makecform('srgb2xyz');
xyzI = applycform(img,xyzT);

%convert xyz to luv colorspace
luvT = makecform('xyz2uvl');
luvI = applycform(xyzI,luvT);


%seperate l,u,v
U=luvI(:,:,1);
L=luvI(:,:,2);
V=luvI(:,:,3);

figure
imshow(luvI, []);
title('LUV colorspace');


%find mean, color variance and color skewness for the channels
extract(1)= mean(U(:));
extract(2) = power(std(U(:)),2);
extract(3)= skewness(U(:));

extract(4) = mean(L(:));
extract(5) = power(std(L(:)),2);
extract(6) = skewness(L(:));

extract(7) = mean(V(:));
extract(8) = power(std(V(:)),2);
extract(9) = skewness(V(:));

extract=extract';
end




function new_hist = edgeDetection(img)
%convert rgb colorspace into ycbcr colorspace to get the luma
img2=rgb2ycbcr(img);


%obtain the Y component
Y=double(img2(:,:,1));

figure
subplot(121)
imshow(img2, []);
title('Ycbcr colorspace');
subplot(122)
imshow(Y, []);
title('Y component');

% sobel filters to be used
T1 = zeros(3,3,5);

T1(:,:,1) = [1 2 1; ...
             0 0 0; ...
            -1 -2 -1]; % Gy (vertical)
        
T1(:,:,2) = [-1 0 1; ...
             -2 0 2; ...
             -1 0 1];   % Gx (horizontal)
         
T1(:,:,3) = [2 2 -1; ...
             2 -1 -1; ...
             -1 -1 -1];% 45 degrees diagonal
         
T1(:,:,4) = [-1 2 2; ...
             -1 -1 2; ...
             -1 -1 -1];% 135 degrees diagonal

T1(:,:,5) = [-1 0 1; ...
              0 0 0; ...
              1 0 -1]; % random

          
          
% apply sobel mask for all T1
for i = 1:5
    sob(:,:,i) = filter2(T1(:,:,i),Y);
end

% calculate the max sobel gradient
[m, n] = max(sob,[],3);

%use canny edge detection to detect edges
edge_d = edge(Y, 'canny');


%edge image multiplied by Sobel mask
im2 = (n.*edge_d);

figure
imshow(im2);
title('Detected edges')

%find hisogram
figure
histogram(im2(:),5)';
title('Image histogram');

new_hist=hist(im2(:),5)';

end
