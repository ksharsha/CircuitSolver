
clc;
clear all;
%%find min i,min,maxi,maxj in each region and label the region accordingly
%%%%%%%%%%%%%%%%%%%%%%%always remember in matlab a(i,j) means i is along y
%%%%%%%%%%%%%%%%%%%%%%%axis and j is along x axis

%%for a single resistor case 0.7*max is working as threshold if not even
%%inductors are getting detected as straight lines

%%the lines appear either for a resistor or for an inductor that's all!!
%%applying connected components to im image

I = imread('download.jpg');
I_plot = I;
if(size(I,3)==3)
    I = rgb2gray(I);
end

%figure ,imshow(I);
%figure('name','Histogram'),imhist(I);
level = graythresh(I);
J = im2bw(I,level);
I = J;
%figure('name','Otsus segmentation'),imshow(J);
dimJ = size(J);
R = zeros(dimJ);
k = 1;
scount = 1;
J(1,:) = 1;
J(:,1) = 1;
for i = 2:dimJ(1)
    for j = 2:dimJ(2)
        if J(i,j) == 0
            if J(i,j-1) == 1 && J(i-1,j) == 1
                R(i,j) = k;
                k = k + 1;
            elseif J(i,j-1) == 1 && J(i-1,j) == 0
                R(i,j) = R(i-1,j);
            elseif J(i,j-1) == 0 && J(i-1,j) == 1
                R(i,j) = R(i,j-1);
            elseif J(i,j-1) == 0 && J(i-1,j) == 0
                R(i,j) = R(i-1,j);
                if R(i,j-1) ~= R(i-1,j)
                    for l = 1:dimJ(1)
                        for m = 1:dimJ(2)
                            if R(l,m) == R(i,j-1)
                                R(l,m) = R(i-1,j);%tracing back
                            end
                        end
                    end
                end
            end
        end
    end
end
%R stores region values
%figure('name','segmented'),imshow(J);%1
sign = [1,1,0];
flag = 0;
for i = 2:dimJ(1)
    for j = 2:dimJ(1)
        if R(i,j) ~= 0
            dimS = size(sign);
            for k = 1:dimS(1)%checks which region the current pixel belongs to and increments that particular value
                if R(i,j) == sign(k,2)
                    sign(k,3) = sign(k,3) + 1;%gives the total number of pixels in that region
                    flag = 1;
                    break
                end
            end
            if flag == 0 %new region
                sign = [sign;dimS(1)+1,R(i,j),1];%adding a new roww to the matrix
            end
            flag = 0;
        end
    end
end
%identifying the region with maximum pixels%
dimS = size(sign);
mx = max(sign(:,3));%%max pixels region

for i = 1:dimS(1)
    if sign(i,3) == mx
        obj = sign(i,2);%%gives the region number
        break
    end
end
%%modifying the program
%%making the mx go to zero
for i = 1:dimS(1)
    if sign(i,3) == mx
        sign(i,3)=1;
        break
    end
end
mx1 = max(sign(:,3));%%second maxima
for i = 1:dimS(1)
    if sign(i,3) == mx1
        obj1 = sign(i,2);
        break
    end
end
if mx1<(mx/2)
    obj1=obj;
end
for i = 1:dimS(1)
    if sign(i,3) == mx1
        sign(i,3)=1;
        break
    end
end
mx2 = max(sign(:,3));%%second maxima
S = zeros(dimJ);
for i = 1:dimJ(1)
    for j = 1:dimJ(2)
        if R(i,j) == obj||R(i,j)== obj1
            S(i,j) = 0;%%making max pixel region to 0 value and rest 1
        else
            S(i,j) = 1;
        end
    end
end
%figure('name','Object Plot'),imshow(S);%2

%%%%%%%%%%%identifying voltage sources%%%%%%%%
[centers,radii,metric]=imfindcircles(S,[10,40],...
    'Sensitivity',0.85)
[center_x,center_y]=size(centers);
flag=zeros(center_x);
%print('the number of identified circles are',center_x);
%%the following for loop checks for false positives
for i=1:center_x
    for j=1:center_x
        if j~=i && pdist2(centers(i,:),centers(j,:))<(radii(i)+radii(j))
            flag(i)=1;
        end
    end
    flag(i)
end
%%if flag variable is set to one it means that that particular circle
%%was a false positive
%%removing these identified voltage source from the image for further
%%processing
a=0;
circle_removed=zeros(dimJ);
circle_removed=S;%%this matrix has the voltage source removed from it
for a=1:center_x
    for i=3:dimJ(1)-2
        for j=3:dimJ(2)-2
            if ((radii(a)-3)<pdist2(centers(a,:),[j,i]))&&(pdist2(centers(a,:),[j,i])<(radii(a)+3))
                circle_removed(i,j)=1-flag(a)*(~S(i,j));%%multiplying with flag to make sure false positives are removed
            end
        end
    end
end
    
   % c2 = bwmorph(circle_removed,'remove');        
        
%[x1, y1] = ginput(1)
%viscircles(centers, radii,'EdgeColor','b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BW = edge(S,'log');
%figure('name','Object Plot'),imshow(circle_removed);hold on;
viscircles(centers, radii,'EdgeColor','red');

m=zeros(dimJ);
for i = 3:dimJ(1)-2
    for j = 3:dimJ(2)-2
       m(i,j)=min(max((abs(circle_removed(i+1,j)-circle_removed(i,j))),(abs(circle_removed(i,j)-circle_removed(i-1,j)))),max((abs(circle_removed(i,j+1)-circle_removed(i,j))),(abs(circle_removed(i,j)-circle_removed(i,j-1)))));
       %min(min((abs(S(i+1,j)-S(i,j))),(abs(S(i,j)-S(i-1,j)))),min((abs(S(i,j+1)-S(i,j))),(abs(S(i,j)-S(i,j-1)))));%%x gradient
       %max((abs(S(i,j+1)-S(i,j))),(abs(S(i,j)-S(i,j-1))));
    end
end
% figure('name','Object Plot'),imshow(m);%5
connected=zeros(dimJ);
a=[1 1 1;1 1 1;1 1 1];
im=imdilate(m,a);
% figure,imshow(im);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%dividing im into regions for localization of components%%
J1=not(im);
dimJ1=size(J1);
im_connected = zeros(dimJ1);
k = 1;
scount = 1;
J1(1,:) = 1;
J1(:,1) = 1;
for i = 2:dimJ1(1)
    for j = 2:dimJ1(2)
        if J1(i,j) == 0
            if J1(i,j-1) == 1 && J1(i-1,j) == 1
                im_connected(i,j) = k;
                k = k + 1;
            elseif J1(i,j-1) == 1 && J1(i-1,j) == 0
                im_connected(i,j) = im_connected(i-1,j);
            elseif J1(i,j-1) == 0 && J1(i-1,j) == 1
                im_connected(i,j) = im_connected(i,j-1);
            elseif J1(i,j-1) == 0 && J1(i-1,j) == 0
                im_connected(i,j) = im_connected(i-1,j);
                if im_connected(i,j-1) ~= im_connected(i-1,j)
                    for l = 1:dimJ1(1)
                        for m1 = 1:dimJ1(2)
                            if im_connected(l,m1) == im_connected(i,j-1)
                                im_connected(l,m1) = im_connected(i-1,j);%tracing back
                            end
                        end
                    end
                end
            end
        end
    end
end
%figure('name','segmented'),imshow(J);%1

sign1=zeros(1,max(im_connected(:)));
for i = 2:dimJ1(1)
    for j = 2:dimJ1(2)
        if im_connected(i,j) ~= 0
           
           sign1(im_connected(i,j))= sign1(im_connected(i,j))+1;
           
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%try to fit a polynomial
BW = edge(circle_removed,'prewitt');
imwrite(BW,'test2.png');
[H,T,R] = hough(BW, 'Theta', 10:0.5:70);  %%% variable1
% figure, imshow(H,[],'XData',T,'YData',R,...
%             'InitialMagnification','fit');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
P  = houghpeaks(H,1000000,'threshold',0.3*max(H(:))); %%% variable2
x = T(P(:,2)); y = R(P(:,1));
% plot(x,y,'s','color','white');
% Find lines and plot them
lines = houghlines(BW,T,R,P,'FillGap',2 ,'MinLength',6); %%% variable3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H1,T1,R1] = hough(BW, 'Theta', -10:-0.5:-70);  %%% variable1
% figure, imshow(H1,[],'XData',T1,'YData',R1,...
%             'InitialMagnification','fit');
% xlabel('\theta1'), ylabel('\rho1');
% axis on, axis normal, hold on;
P1  = houghpeaks(H1,1000000,'threshold',0.3*max(H1(:))); %%% variable2
x1 = T1(P1(:,2)); y1 = R1(P1(:,1));
% plot(x1,y1,'s','color','white');
% Find lines and plot them
lines1 = houghlines(BW,T1,R1,P1,'FillGap',2 ,'MinLength',6); %%% variable3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure, imshow(S), hold on
count=0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
    if (im(xy(1,2),xy(1,1))==1 && im(xy(2,2),xy(2,1))==1)%%removes false positives
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','blue');
    plot(xy(2,1),xy(2,2),'o','LineWidth',2,'Color','red');
    end
   
  
end
hold on;
region=[1000];
count1=0;
flag3=0;
for k = 1:length(lines1)
   xy = [lines1(k).point1; lines1(k).point2];
    if (im(xy(1,2),xy(1,1))==1 && im(xy(2,2),xy(2,1))==1)%%removes false positives
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','blue');
    plot(xy(2,1),xy(2,2),'o','LineWidth',2,'Color','red');
    end
    if ((im_connected(xy(1,2),xy(1,1))==im_connected(xy(2,2),xy(2,1)))&&(im_connected(xy(1,2),xy(1,1))~=0))%%making sure the line belongs to the same region
        for count1=1:length(region)
            if(im_connected(xy(1,2),xy(1,1))==region(count1))
                flag3=1;
            end
        end
        if(flag3==0)
            region=[region;im_connected(xy(1,2),xy(1,1))];
        end
        flag3=0;
       % region=[region;im_connected(xy(1,2),xy(1,1))];
        
        
    end
end
resistors(1:length(region)-1)=region(2:length(region));
length_res=length(resistors);
pixel_resistors=resistors;
index_resistors=resistors;

i=1;
j=1;
[s1,s2]=size(sign1);
for i=1:length_res     
    pixel_resistors(i)=sign1(resistors(i));  
end
final1pixel_resistors=[1000];
index1pixel_resistors=[1000];
mz=max(pixel_resistors);
for i=1:length_res
    if pixel_resistors(i)>(mz/1.5)%%assuming every resistor is atleast 1/3 the max resistor
        final1pixel_resistors=[final1pixel_resistors;pixel_resistors(i)];
        index1pixel_resistors=[index1pixel_resistors;resistors(i)];
        
    end
end
pixel_res=length(final1pixel_resistors);
finalpixel_resistors(1:pixel_res-1)=final1pixel_resistors(2:pixel_res);
indexpixel_resistors(1:pixel_res-1)=index1pixel_resistors(2:pixel_res);
        
display('The number of resistors idenitifed in the circuit are');
display(length(finalpixel_resistors));
display('The indices of the idenitifed resistors in the circuit are');
display(indexpixel_resistors);
display('The number of pixels in each index are');
display(finalpixel_resistors);
hold on;

matrix_index=zeros(length(finalpixel_resistors),5);
for i=1:length(finalpixel_resistors)
    matrix_index(i,1)=indexpixel_resistors(i);
    matrix_index(i,2)=0;
    matrix_index(i,3)=0;
    matrix_index(i,4)=1000;
    matrix_index(i,5)=1000;
end
c_removed=circle_removed;
%%matrix_index(a,b,c,d,e) a=region,b=maxy,c=maxx,d=miny,e=minx
for i=1:dimJ1(1)
    for j=1:dimJ(2)
        if im(i,j)==1
            for k=1:length(finalpixel_resistors)
                if im_connected(i,j)==matrix_index(k,1)
                    if matrix_index(k,2)<i
                        matrix_index(k,2)=i;
                    end
                    if matrix_index(k,3)<j
                        matrix_index(k,3)=j;
                    end
                    if matrix_index(k,4)>i
                        matrix_index(k,4)=i;
                    end
                    if matrix_index(k,5)>j
                        matrix_index(k,5)=j;
                    end
                    c_removed(i,j)=1;
                end
            end
        end
    end
end
matrix_index % 1st column(ignore), next 4 columns contains co-ods of diagonal of resistor window
%figure;imshow(im_connected);
%                     
% imwrite(c_removed,'comp_removed.jpg');
% a = imread('comp_removed.jpg');

I = ident2(I_plot, matrix_index, centers, radii, c_removed);
                
        




            
                






