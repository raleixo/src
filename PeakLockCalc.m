function pk=PeakLockCalc(u,nbins)
%PEAKLOCKCALC Auxiliary function to PeakLock. It calculates the histograms
%and the peak lock value. It plots all histograms too.
%   Calculates histogams of the velocity matrix.

    %[u,nbins]=deal(vargin{:});

    xi=size(u,1); yi=size(u,2);
    u11=reshape(u,xi*yi,1);
    u11=u11(isfinite(u11)); %filtering nan
    
    [n1,b1]=hist(u11,nbins);
    
    %Vmod1
    u11=abs(u11);
    u11frac1=u11-floor(u11);
    [n2,b2]=hist(u11frac1,nbins);
    
    %Vomd05
    a1=find(u11frac1>0.5);
    u11frac2=1-u11frac1(a1);
    [n3,b3]=hist(u11frac2,nbins);
    
%     figure, bar(b3,n3/sum(n3));
    CenterOfMass=sum((n3.*b3)/sum(n3));
    
    PKLock=4*(0.25-CenterOfMass);
    titlefigure=['Peak Lock = ' num2str(PKLock)];
    figure, subplot(2,2,1), bar(b1,n1/sum(n1)), hold on,...
            subplot(2,2,2), bar(b2,n2/sum(n2)), hold on,...
            subplot(2,2,3), bar(b3,n3/sum(n3)); xlim([0 0.5]); title(titlefigure)
        hold off
        
 pk=PKLock;
        
    
