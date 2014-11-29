function [x0,y0]=intpeak(x1,y1,R,method,N)
result_conv=R;
a=rand(1);
save(['C:\RAleixo\Imagens\Pure' '\' 'RRR' num2str(a) '.mat'],'R','x1','y1','N');
clear R
R=result_conv(y1,x1);
Rxm1=result_conv(y1,x1-1); Rxp1=result_conv(y1,x1+1);
Rym1=result_conv(y1-1,x1); Ryp1=result_conv(y1+1,x1);

% INTPEAK - interpolate correlation peaks in PIV
%
% TimeStamp: November 5th 2014
% Changed by AMRicardo
%
% function [x0,y0]=intpeak(x1,x2,x3,y1,y2,y3,method,N)
% METHOD = 
% 1 for centroid fit, 
% 2 for gaussian fit, 
% 3 for parabolic fit,
% 4 for 2D gaussian fit %Added by AMRicardo
% x1 and y1 are maximal values in respective directions.
% N is interrogation window size. N is either 1x1 or 1x2
%
% This is a subfunction to MATPIV

% Time stamp: 12:32, Apr. 14, 2004.
%
% Copyright 1998-2004, J. Kristian Sveen, 
% jks@math.uio.no/jks36@damtp.cam.ac.uk
% Dept of Mathmatics, University of Oslo/ 
% DAMTP, Univ. of Cambridge, UK
%
% Distributed under the GNU general public license.
%
% For use with MatPIV 1.6 and subsequent versions

if length(N)==2
    M=N(1); N=N(2);
else
    M=N;
end

if any(find(([R Rxm1 Rxp1 Rym1 Ryp1])==0))
    % to avoid Log of Zero warnings
    method=1;
end

if method==1  
    x01=(((x1-1)*Rxm1)+(x1*R)+((x1+1)*Rxp1)) / (Rxm1+ R+Rxp1);
    y01=(((y1-1)*Rym1)+(y1*R)+((y1+1)*Ryp1)) / (Rym1+ R+Ryp1);
    x0=x01-(M);               
    y0=y01-(N);
elseif method==2  
    x01=x1 + ( (log(Rxm1)-log(Rxp1))/( (2*log(Rxm1))-(4*log(R))+(2*log(Rxp1))) );
    y01=y1 + ( (log(Rym1)-log(Ryp1))/( (2*log(Rym1))-(4*log(R))+(2*log(Ryp1))) );  
    x0=x01-(M);
    y0=y01-(N);  
elseif method==3
    x01=x1 + ( (Rxm1-Rxp1)/( (2*Rxm1)-(4*R)+(2*Rxp1)) );
    y01=y1 + ( (Rym1-Ryp1)/( (2*Rym1)-(4*R)+(2*Ryp1)) ); 
    x0=x01-(M);
    y0=y01-(N);
    
elseif method==4
    y=y1; x=x1;
    c10=zeros(3,3);
    c01=c10;c11=c10;c20=c10;c02=c10;
     for i=-1:1
        for j=-1:1
            %following 15 lines based on
            %H. Nobach Æ M. Honkanen (2005)
            %Two-dimensional Gaussian regression for sub-pixel displacement
            %estimation in particle image velocimetry or particle position
            %estimation in particle tracking velocimetry
            %Experiments in Fluids (2005) 38: 511–515
            c10(j+2,i+2)=i*log(result_conv(y+j, x+i));
            c01(j+2,i+2)=j*log(result_conv(y+j, x+i));
            c11(j+2,i+2)=i*j*log(result_conv(y+j, x+i));
            c20(j+2,i+2)=(3*i^2-2)*log(result_conv(y+j, x+i));
            c02(j+2,i+2)=(3*j^2-2)*log(result_conv(y+j, x+i));
            %c00(j+2,i+2)=(5-3*i^2-3*j^2)*log(result_conv_norm(maxY+j, maxX+i));
        end
    end
    c10=(1/6)*sum(sum(c10));
    c01=(1/6)*sum(sum(c01));
    c11=(1/4)*sum(sum(c11));
    c20=(1/6)*sum(sum(c20));
    c02=(1/6)*sum(sum(c02));
    %c00=(1/9)*sum(sum(c00));
    
    deltax=(c11*c01-2*c10*c02)/(4*c20*c02-c11^2);
    deltay=(c11*c10-2*c01*c20)/(4*c20*c02-c11^2);
    x01=x+deltax;
    y01=y+deltay;
    
% %     SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
% %     SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
% %     vector=[SubpixelX, SubpixelY];
    x0=x01-(M);
    y0=y01-(N);
    
elseif method==5
    [Xo,Yo]=meshgrid(-2:1:2);
    [Xq,Yq]=meshgrid(-1:0.01:1);
    Vo=result_conv(x1-2:1:x1+2,y1-2:1:y1+2);
    Vq = interp2(Xo,Yo,Vo,Xq,Yq,'cubic');
    x01=x1+Xq(Vq==max(max(Vq)));
    y01=y1+Yq(Vq==max(max(Vq)));
    if length(x01)~=1, x01=0; y01=0; end
    x0=x01-(M);
    y0=y01-(N);
else
    
    disp(['Please include your desired peakfitting function; 1 for',...
	  ' 3-point fit, 2 for gaussian fit, 3 for parabolic fit'])
    
end

x0=real(x0);
y0=real(y0);


%% PIVLAB
%(result_conv,interrogationarea,x,y,SubPixOffset);
% result_conv: normalized (max=255) correlation matrice, has the size of the interrogation area
% result_conv =fftshift(real(ifft2(conj(fft2(image1_crop)).*fft2(image2_crop))));
% result_conv=result_conv/max(max(result_conv))*255; %normalize, peak=always 255
% [y,x] = find(result_conv==255); %Find the 255 peak (it takes the first
% peak if there is more than one
% SubPixOffset: Remainder after division (rem)
% if (rem(interrogationarea,2) == 0) %for the subpixel displacement measurement
%     SubPixOffset=1;
% else
%     SubPixOffset=0.5;
% end

% SUBPIXGAUSS: 3x2 point gaussian (the same implemented in MatPiv

% function [vector] = SUBPIXGAUSS (result_conv,interrogationarea,x,y,SubPixOffset)
% if (x <= (size(result_conv,1)-1)) && (y <= (size(result_conv,1)-1)) && (x >= 1) && (y >= 1)
%     %the following 8 lines are copyright (c) 1998, Uri Shavit, Roi Gurka, Alex Liberzon, Technion – Israel Institute of Technology
%     %http://urapiv.wordpress.com
%     f0 = log(result_conv(y,x));
%     f1 = log(result_conv(y-1,x));
%     f2 = log(result_conv(y+1,x));
%     peaky = y+ (f1-f2)/(2*f1-4*f0+2*f2);
%     f0 = log(result_conv(y,x));
%     f1 = log(result_conv(y,x-1));
%     f2 = log(result_conv(y,x+1));
%     peakx = x+ (f1-f2)/(2*f1-4*f0+2*f2);
%     %
%     SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
%     SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
%     vector=[SubpixelX, SubpixelY];  
%  else
%     vector=[NaN NaN];
% end
% function [vector] = SUBPIX2DGAUSS (result_conv,interrogationarea,x,y,SubPixOffset)
% if (x <= (size(result_conv,1)-1)) && (y <= (size(result_conv,1)-1)) && (x >= 1) && (y >= 1)
%     c10=zeros(3,3);
%     c01=c10;c11=c10;c20=c10;c02=c10;
%     for i=-1:1
%         for j=-1:1
%             %following 15 lines based on
%             %H. Nobach Æ M. Honkanen (2005)
%             %Two-dimensional Gaussian regression for sub-pixel displacement
%             %estimation in particle image velocimetry or particle position
%             %estimation in particle tracking velocimetry
%             %Experiments in Fluids (2005) 38: 511–515
%             c10(j+2,i+2)=i*log(result_conv(y+j, x+i));
%             c01(j+2,i+2)=j*log(result_conv(y+j, x+i));
%             c11(j+2,i+2)=i*j*log(result_conv(y+j, x+i));
%             c20(j+2,i+2)=(3*i^2-2)*log(result_conv(y+j, x+i));
%             c02(j+2,i+2)=(3*j^2-2)*log(result_conv(y+j, x+i));
%             %c00(j+2,i+2)=(5-3*i^2-3*j^2)*log(result_conv_norm(maxY+j, maxX+i));
%         end
%     end
%     c10=(1/6)*sum(sum(c10));
%     c01=(1/6)*sum(sum(c01));
%     c11=(1/4)*sum(sum(c11));
%     c20=(1/6)*sum(sum(c20));
%     c02=(1/6)*sum(sum(c02));
%     %c00=(1/9)*sum(sum(c00));
%     
%     deltax=(c11*c01-2*c10*c02)/(4*c20*c02-c11^2);
%     deltay=(c11*c10-2*c01*c20)/(4*c20*c02-c11^2);
%     peakx=x+deltax;
%     peaky=y+deltay;
%     
%     SubpixelX=peakx-(interrogationarea/2)-SubPixOffset;
%     SubpixelY=peaky-(interrogationarea/2)-SubPixOffset;
%     vector=[SubpixelX, SubpixelY];
% else
%     vector=[NaN NaN];
% end