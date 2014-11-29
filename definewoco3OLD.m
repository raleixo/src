function [comap,A1,world]=definewoco3(filename,typ,path)
% DEFINEWOCO - calculate the mapping from image to physical coordinates in MatPIV
%
% [comap,A1,world]=definewoco(image,coordinatestyle, path)
% 
% DEFINEWOCO2 is a file for calculating the pixel to world 
% coordinate transformation. Calibration file is saved in path

% This function needs an image with distinct coordinate points. 
% These points should either be easily locatable dots in the 
% image (these are assumed to be circular), or they will have to
% be defined by a grid with known spacing. Your points need to be
% WHITE on a BLACK (dark) background.
%
% definewoco2('woco.bmp','o'); will assume circular points with a well 
% defined peak
%
% definewoco2('woco.bmp','+'); will assume that the image contains a grid 
% with known spacing. In this case the image is cross-correlated with an 
% artificial cross.
%
% definewoco2('woco.bmp','x'); will assume that the image contains grid 
% points in the form of x's. In this case the image is cross-correlated 
% with an artificial x (rotated cross).
%
% definewoco2('woco.bmp','.'); will assume circular points where the peak 
% is not well defined, for example if it is flat and relatively large. In 
% this case the input-image is cross correlated with a gaussian bell to 
% emphasize the center of each point in the original image. This option has 
% now the added functionality of letting the user enter an aproximate size
% of his/her world-coordinate point. This helps in many cases where the 
% points are "very" wide, for example 20 pixels in diameter.
%
% definewoco2('woco.bmp','s'); will assume that the image contains a checkboard
%In this case the image is cross-correlated 
% with an artificial _| .
%
% The user will then have to mark the local regions around each 
% coordinate point using the mouse (first button).
%
% Subsequently the user will be asked to enter the physical coordinates
% (cm) for each point.
%
% In the final step one will have to choose between a linear and non-
% linear mapping function to use in the calculation of the mapping 
% factors. In most cases a linear function will do a proper job.
%
% The final result will be saved to a file in the present work directory. 
% The file is named 'worldcoX.mat', where the X is any number specified 
% by the user. This option is for the cases where one might have two or 
% more different coordinate systems in the same working directory. If this 
% is not the case just press <enter> when asked for a number. The file will 
% then be named 'worldco.mat'
%
% See also: MATPIV


% Copyright, J. Kristian Sveen, 1999-2004, last revision April 16, 2004
%             jks@math.uio.no
% For use with MatPIV 1.6.1
% Distributed under the GNU general public license
format long
%filename=[seqcalib.path seqcalib.name];
typ='.';
if ischar(filename)
  A=imread(filename);
else
  A=filename;
end

if size(A,3)>1, A=double(rgb2gray(A)); else A=double(A); end

A=255-A;
[ay,ax]=size(A);

my_ver=version;
my_ver=str2num(my_ver(1:3));
% if my_ver>=6.5, pixval on, end

if strcmp(typ,'+')==1
  load articross.mat
  disp('....calculating....this may take a few seconds.')
  b=A./max(A(:)); % normalize in order to make max(A(:))=1 (since max(kr(:))=1);
  A=xcorrf2(b-mean(b(:)),kr-mean(kr(:)))./(std(kr(:))*std(b(:))*size(kr,1)*ay);
  [ax,ay]=size(A); [bx,by]=size(b);
  dx=(ax-bx+1)/2; dy=(ay-by+1)/2;
  A=A(dy+1:end-(dy-1),dx+1:end-(dx-1));
  disp('...Done!')
  disp('Now mark the crosses you whish to use as coordinate points')
elseif strcmp(typ,'x')==1
  load articross2.mat
  kr=double(kr);
  disp('....calculating....this may take a few seconds.')
  b=A./max(A(:)); % normalize in order to make max(A(:))=1 (since max(kr(:))=1);
  A=xcorrf2(b-mean(b(:)),kr-mean(kr(:)))./(std(kr(:))*std(b(:))*size(kr,1)*ay);
  [ax,ay]=size(A); [bx,by]=size(b);
  dx=(ax-bx+1)/2; dy=(ay-by+1)/2;
  A=A(dy+1:end-(dy-1),dx+1:end-(dx-1));
  disp('...Done!')
  disp('Now mark the crosses you whish to use as coordinate points')
elseif strcmp(typ,'o')==1
  %no need to do anything with the image in this case
elseif strcmp(typ,'.')==1
%   disp('Please give the approximate width of your points (in pixels -')
%   point_size=input(['default is 20). Type 0 here to get old behaviour of definewoco: ']);
%   if isempty(point_size), point_size=30; 
%   else, point_size=point_size+10; end
 point_size=100+10;
  w=weight('gaussian',point_size,0.1);
  disp('....calculating....this may take a few seconds.')
  b=A./max(A(:)); % normalize in order to make max(A(:))=1 (since max(w(:))=1);
  A=xcorrf2(b-mean(b(:)),w-mean(w(:)))./...
    (std(b(:))*std(w(:))*(point_size+10)*ay);
  A=A(point_size/2 +1:end-(point_size/2 -1),...
      point_size/2 +1:end-(point_size/2 -1));
  disp('...Done!')
  disp('Now mark the dots you whish to use as coordinate points')
  elseif strcmp(typ,'x')==1
  load articross2.mat
  kr=double(kr);
  disp('....calculating....this may take a few seconds.')
  b=A./max(A(:)); % normalize in order to make max(A(:))=1 (since max(kr(:))=1);
  A=xcorrf2(b-mean(b(:)),kr-mean(kr(:)))./(std(kr(:))*std(b(:))*size(kr,1)*ay);
  [ax,ay]=size(A); [bx,by]=size(b);
  dx=(ax-bx+1)/2; dy=(ay-by+1)/2;
  A=A(dy+1:end-(dy-1),dx+1:end-(dx-1));
  disp('...Done!')
  disp('Now mark the crosses you whish to use as coordinate points')
elseif strcmp(typ,'s')==1
  load artisquare4.mat
  kr=double(kr);
  disp('....calculating....this may take a few seconds.')
  b=A./max(A(:)); % normalize in order to make max(A(:))=1 (since max(kr(:))=1);
  A=xcorrf2(b-mean(b(:)),kr-mean(kr(:)))./(std(kr(:))*std(b(:))*size(kr,1)*ay);
  [ax,ay]=size(A); [bx,by]=size(b);
  dx=(ax-bx+1)/2; dy=(ay-by+1)/2;
  A=A(dy+1:end-(dy-1),dx+1:end-(dx-1));
  disp('...Done!')
  disp('Now mark the crosses you whish to use as coordinate points')

else
  disp('Not a valid coordinate style. Please use either a + or an o')
  return
end
%save('A11.mat','A')
figure
imagesc(A), axis('equal');
usr1=1;

disp('Please mark your world coordinate points with left mouse button.');
disp('Press ENTER when finished!')

[x1,y1]=mginput;
x1=round(x1); y1=round(y1);
for i=1:1:size(x1,1)
    if y1(i)-9<1, edgy1=(y1(i)-9)-1; else edgy1=0; end
    if y1(i)+8>size(A,1), edgy2=(y1(i)+8)-size(A,1); else edgy2=0; end
    
    if x1(i)-9<1, edgx1=(x1(i)-9)-1; else edgx1=0;  end
    if x1(i)+8>size(A,2), edgx2=(x1(i)+8)-size(A,2); else edgx2=0; end
    
    B=A( y1(i)-9+edgy1:y1(i)+8+edgy2, x1(i)-9+edgx1:x1(i)+8+edgx2);
    %figure(2), imagesc(B)
    [max_y1 max_x1]=find(B==max(max(B)));
    if size(max_x1,1)>1 
        max_x1=max_x1(2);   max_y1=max_y1(2);
    end  
    [x0 y0]=intpeak(max_x1,max_y1,B(max_y1,max_x1),...
        B(max_y1,max_x1-1),B(max_y1,max_x1+1),...
        B(max_y1-1,max_x1),B(max_y1+1,max_x1),1,9);
    x(i)=x1(i)+x0;
    y(i)=y1(i)+y0;
end
newx=x(3);
nx=size(A,2); ny=size(A,1);
dx0=x(3)-x(1); dy0=abs(y(2)-y(1));
dx=dx0;
p=1;
clear x1 y1;
newx1=newx;

%this is for line from x(3) up to end of image
while newx<(nx-2*dx0)
   
    newxaux=newx+dx;
    Naux(i)=newxaux;
    newy=y(3);
    %check if point is over dot
    y1=round(newy); x1=round(newxaux);
    i=1;
    if y1-9<1, edgy1=(y1-9)-1; else edgy1=0; end
    if y1+8>size(A,1), edgy2=(y1+8)-size(A,1); else edgy2=0; end
    
    if x1-9<1, edgx1=(x1-9)-1; else edgx1=0;  end
    if x1+8>size(A,2), edgx2=(x1+8)-size(A,2); else edgx2=0; end
    
    B=A( y1-9+edgy1:y1+8+edgy2, x1-9+edgx1:x1+8+edgx2);
    %figure(2), imagesc(B)
    [max_y1 max_x1]=find(B==max(max(B)));
    if size(max_x1,1)>1 
        max_x1=max_x1(2);   max_y1=max_y1(2);
    end  
    
    if (max_x1<17 && max_y1<17)
       
        [x0 y0]=intpeak(max_x1,max_y1,B(max_y1,max_x1),...
            B(max_y1,max_x1-1),B(max_y1,max_x1+1),...
            B(max_y1-1,max_x1),B(max_y1+1,max_x1),1,9);
        x2(p)=x1+x0;
        y2(p)=y1+y0;
        dx=x2(p)-newx;
        newx=(x2(p));
        
        p=p+1;
    else
        dx=dx+0.1*dx;
    end
%     newx1=newx1+0;
end
 
%this is for line from x(3) to begining of image 1
newx=x(3)+dx0;
dx=dx0;
while newx>(1+4*dx0)
   
    newxaux=newx-dx;
    Naux(i)=newxaux;
    newy=y(3);
    %check if point is over dot
    y1=round(newy); x1=round(newxaux);
    i=1;
    if y1-9<1, edgy1=(y1-9)-1; else edgy1=0; end
    if y1+8>size(A,1), edgy2=(y1+8)-size(A,1); else edgy2=0; end
    
    if x1-9<1, edgx1=(x1-9)-1; else edgx1=0;  end
    if x1+8>size(A,2), edgx2=(x1+8)-size(A,2); else edgx2=0; end
    
    B=A( y1-9+edgy1:y1+8+edgy2, x1-9+edgx1:x1+8+edgx2);
    %figure(2), imagesc(B)
    [max_y1 max_x1]=find(B==max(max(B)));
    if size(max_x1,1)>1 
        max_x1=max_x1(2);   max_y1=max_y1(2);
    end  
    
    if (max_x1<17 && max_y1<17)
         
        [x0 y0]=intpeak(max_x1,max_y1,B(max_y1,max_x1),...
            B(max_y1,max_x1-1),B(max_y1,max_x1+1),...
            B(max_y1-1,max_x1),B(max_y1+1,max_x1),1,9);
        x3(p)=x1+x0;
        y3(p)=y1+y0;
        dx=abs(x3(p)-newx);
        newx=(x3(p));
        
        p=p+1;
    else
        dx=dx+0.1*dx;
    end
%     newx1=newx1+0;
end

%this is for column from y(2) to the bottom of the image
dy=dy0;
newy=y(2);
while newy<(size(A,1)-4*dy0)
   
    newyaux=newy+dy;
    Nauxy(i)=newyaux;
    newx=x(2);
    %check if point is over dot
    x1=round(newx); y1=round(newyaux);
    i=1;
    if y1-9<1, edgy1=(y1-9)-1; else edgy1=0; end
    if y1+8>size(A,1), edgy2=(y1+8)-size(A,1); else edgy2=0; end
    
    if x1-9<1, edgx1=(x1-9)-1; else edgx1=0;  end
    if x1+8>size(A,2), edgx2=(x1+8)-size(A,2); else edgx2=0; end
    
    B=A( y1-9+edgy1:y1+8+edgy2, x1-9+edgx1:x1+8+edgx2);
    %figure(2), imagesc(B)
    [max_y1 max_x1]=find(B==max(max(B)));
    if size(max_x1,1)>1 
        max_x1=max_x1(2);   max_y1=max_y1(2);
    end  
    
    if (max_x1<17 && max_y1<17)
       
        [x0 y0]=intpeak(max_x1,max_y1,B(max_y1,max_x1),...
            B(max_y1,max_x1-1),B(max_y1,max_x1+1),...
            B(max_y1-1,max_x1),B(max_y1+1,max_x1),1,9);
        x4(p)=x1+x0;
        y4(p)=y1+y0;
        dx=abs(y4(p)-newy);
        newy=(y4(p));
        
        p=p+1;
    else
        dy=dy+0.1*dy;
    end
%     newx1=newx1+0;
end

dy=dy0;
newy=y(2)+dy0;
while newy>(1+5*dy0)
   
    newyaux=newy-dy;
    Nauxy(i)=newyaux;
    newx=x(2);
    %check if point is over dot
    x1=round(newx); y1=round(newyaux);
    i=1;
    if y1-9<1, edgy1=(y1-9)-1; else edgy1=0; end
    if y1+8>size(A,1), edgy2=(y1+8)-size(A,1); else edgy2=0; end
    
    if x1-9<1, edgx1=(x1-9)-1; else edgx1=0;  end
    if x1+8>size(A,2), edgx2=(x1+8)-size(A,2); else edgx2=0; end
    
    B=A( y1-9+edgy1:y1+8+edgy2, x1-9+edgx1:x1+8+edgx2);
    %figure(2), imagesc(B)
    [max_y1 max_x1]=find(B==max(max(B)));
    if size(max_x1,1)>1 
        max_x1=max_x1(2);   max_y1=max_y1(2);
    end  
    
    if (max_x1<17 && max_y1<17)
         
        [x0 y0]=intpeak(max_x1,max_y1,B(max_y1,max_x1),...
            B(max_y1,max_x1-1),B(max_y1,max_x1+1),...
            B(max_y1-1,max_x1),B(max_y1+1,max_x1),1,9);
        x5(p)=x1+x0;
        y5(p)=y1+y0;
        dx=abs(y5(p)-newy);
        newy=(y5(p));
        
        p=p+1;
    else
        dy=dy+0.1*dy;
    end
%     newx1=newx1+0;
end

XG=[x2 x3]; YG=[y4 y5];
XG=sort(XG(find(XG)),'ascend')';
YG=sort(YG(find(YG)),'ascend')';
[aa bb]=meshgrid(XG,YG);

figure, imagesc(A), hold on, plot(aa,bb,'k+'), axis equal
%save(['G:\WeimingWuPC\RAleixo\' 'coords.mat'], 'x', 'y', 'nx', 'ny', 'x2', 'y2', 'Naux', 'dx', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5');

%close(2)
disp('Now you need to give the physical coordinates to each of the points specified!')
disp('-----------------------')
hold on

worldx=((1:size(XG,1))-1)*0.01;
worldy=((1:size(YG,1))-1)*0.01;
[wcoordx wcoordy]=meshgrid(worldx,worldy);
aa0=aa(1);
aa=aa-aa0;
bb0=bb(1);
bb=bb-bb0;

p=1;
dxtarget=0.01; dytarget=0.01; %coordinates 
for i=1:size(aa,2)
    for j=1:size(aa,1)
        [world1(p,1:2)]=[(i-1)*dxtarget (j-1)*dytarget];
        [wcc(p,1:2)]=[aa(1,i) bb(j,1)];
        p=p+1;
    end
end
            
% for i=1:1:size(wcc,1)
%   hs=plot(wcc(i,1)+aa0,wcc(i,2)+bb0,'wo');
%   set(hs,'MarkerSize',[16])
%   [world(i,1:2)]=input(['Please enter the world coordinates for the white \n circle  marked in the current figure (in square parenthesis): ']);
%   set(hs,'MarkerFaceColor',[0.1 0.1 0.1])
% end
% x=world1(:,1); y=world1(:,2);
% world=wcc;
x=wcc(:,1); y=wcc(:,2);
world=world1;


% Construct function for fitting.
inpx='inne';
while strcmp(inpx,'inne')
  mapfun=input('Mapping function. (N)onlinear or (L)inear (N/n/L/l): ','s');
  if strcmp(mapfun,'N')==1 || strcmp(mapfun,'n')==1
    if length(world)>=6
       A1=[ones(size(x,1),1) x y (x.*y) (x.^2) (y.^2)];%A1=[ones(size(x,1),1) x.' y.' (x.*y).' (x.^2).' (y.^2).'];
      inpx='ute';
    else
      disp('Not enough points specified to calculate nonlinear mapping factors.')
      disp('Using linear mapping function.');
      A1=[ones(size(x,1),1) x y]; 
      inpx='ute';
    end
  elseif strcmp(mapfun,'L')==1 || strcmp(mapfun,'l')==1
    A1=[ones(size(x,1),1) x y];
    inpx='ute';
  else
    disp('Please specify mapping function! (N/n/L/l)')
  end
end

comap=(A1\world(:,:));  % Fit using a minimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for error
err=norm( (A1*comap-world));
disp(['Error (norm) = ',num2str(err)])
% give a warning if error is larger than a certain threshold
% 1 is just chosen as a test case. This needs testing.
if err>1
  disp(['WARNING! The minimized system of equations has a large ', ...
	'error.'])
  disp('Consider checking your world coordinate input')
  if strcmp(mapfun,'L') || strcmp(mapfun,'l')
    disp('Also consider using a nonlinear mapping function');
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inp=input('Save coordinate file....specify number >> ');
navn=['worldco',num2str(inp)];
save([path '\' navn],'comap')
disp(['Coordinate mapping factors saved to file:  ', [path '\' navn] ])

close %close window containing the image

