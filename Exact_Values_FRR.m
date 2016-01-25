function [xOpt, yOpt]=Exact_Values_FRR(x, y, Pressure, Time, Frequency)
%function [xOpt, yOpt, Legendtstart]=Exact_Values_FRR(x, y, Pressure, Time, Frequency)

Pressure = round(Pressure,1); %round to tenths to make process less sensitive.  Adjust rounding to adjust sensitivity

b = mod(length(x),2);

if b == 0
    
    xstartIndex = find(Time < x(1),1,'last'); %find first time input value.  The input most likely won't match exactly to time values
    
    %Check starting index catches start of fill
    while Pressure(xstartIndex) > Pressure(xstartIndex - Frequency) %go left on x axis
        xstartIndex = xstartIndex - Frequency;
    end
    
    step = ceil(Frequency*0.5);
    %Confirm with pressure rise
    while (Pressure(xstartIndex+1*step*step) >= Pressure(xstartIndex+2*step)||Pressure(xstartIndex) <= Pressure(xstartIndex + 1*step))...
            && (Pressure(xstartIndex+2*step) >= Pressure(xstartIndex+3*step)||Pressure(xstartIndex) >= Pressure(xstartIndex + 1*step))...
            && (Pressure(xstartIndex+3*step) >= Pressure(xstartIndex+4*step)||Pressure(xstartIndex) >= Pressure(xstartIndex + 1*step))...
            && (Pressure(xstartIndex+4*step) >= Pressure(xstartIndex+5*step)||Pressure(xstartIndex) >= Pressure(xstartIndex + 1*step))...
            && (Pressure(xstartIndex+5*step) >= Pressure(xstartIndex+6*step)||Pressure(xstartIndex) >= Pressure(xstartIndex + 1*step))...
            && (Pressure(xstartIndex+6*step) >= Pressure(xstartIndex+7*step)||Pressure(xstartIndex) >= Pressure(xstartIndex + 1*step))%go right on x axis
        xstartIndex = xstartIndex + 1*step;
    end
    
    xOpt(1) = Time(xstartIndex);
    yOpt(1) = Pressure(xstartIndex);
    
    %% optimize all tLeakBeginning values
    for i=2:2:length(x)-1
        
        %         IndexArray = x(i,1)>Time;                                       %Array of Time Values <x(1,1)
        %         xstartIndex = numel(IndexArray(IndexArray))-1;                  %Number of Time Values <x(1,1)
        xstartIndex = find(Time < x(i),1,'last'); %find first time input value.  The input most likely won't match exactly to time values
        
        %Check starting index start of leak check
        while Pressure(xstartIndex) <= Pressure(xstartIndex - Frequency) %go left on x axis
            xstartIndex = xstartIndex - Frequency;
        end
        while Pressure(xstartIndex) > Pressure(xstartIndex - 3)          %go right on x axis
            xstartIndex = xstartIndex + 1;
        end
        
        xOpt(i) = Time(xstartIndex-2);
        yOpt(i) = Pressure(xstartIndex-2);
    end
    
    %% tLeakEnd
    for i=3:2:length(x)-1
        
        xstartIndex = find(Time < x(i),1,'last'); %find first time input value.  The input most likely won't match exactly to time values
        
        while Pressure(xstartIndex) >= Pressure(xstartIndex + Frequency) %go right on x axis
            xstartIndex = xstartIndex + Frequency;
        end
        while Pressure(xstartIndex) < Pressure(xstartIndex + Frequency)  %go left on x axis
            xstartIndex = xstartIndex - Frequency;
        end
        
        xOpt(i) = Time(xstartIndex + Frequency);
        yOpt(i) = Pressure(xstartIndex + Frequency);
    end
    
    %% tEnd
    xend = length(x);
    xstartIndex = find(Time < x(xend),1,'last');
    
    if Pressure(xstartIndex) <= Pressure(xstartIndex-3)             %go left on x axis
        while Pressure(xstartIndex) <= Pressure(xstartIndex-3)
            xstartIndex = xstartIndex-4*Frequency;
        end
        xOpt(xend) = Time(xstartIndex);
        yOpt(xend) = Pressure(xstartIndex);
    end
    if Pressure(xstartIndex) > Pressure(xstartIndex-2*Frequency)    %go right on x axis
        while Pressure(xstartIndex)<max(Pressure) && xstartIndex<numel(Pressure)
            xstartIndex = xstartIndex+1;
        end
        xOpt(xend) = Time(xstartIndex);
        yOpt(xend) = Pressure(xstartIndex);
    end
    
else
    %% t0
    xstartIndex = find(Time < x(i),1,'last');

    while Pressure(xstartIndex) > Pressure(xstartIndex-Frequency)       %go left on x axis
        xstartIndex = xstartIndex -Frequency;
    end
    while (Pressure(xstartIndex+2) <= Pressure(xstartIndex+1)||Pressure(xstartIndex+1) <= Pressure(xstartIndex))...
            && (Pressure(xstartIndex+3) <= Pressure(xstartIndex+2)||Pressure(xstartIndex+1) <= Pressure(xstartIndex))...
            && (Pressure(xstartIndex+4) <= Pressure(xstartIndex+3)||Pressure(xstartIndex+1) <= Pressure(xstartIndex))...
            && (Pressure(xstartIndex+5) <= Pressure(xstartIndex+4)||Pressure(xstartIndex+1) <= Pressure(xstartIndex))...
            && (Pressure(xstartIndex+6) <= Pressure(xstartIndex+5)||Pressure(xstartIndex+1) <= Pressure(xstartIndex))...
            && (Pressure(xstartIndex+7) <= Pressure(xstartIndex+6)||Pressure(xstartIndex+1) <= Pressure(xstartIndex))%go right on x axis
        xstartIndex = xstartIndex + 1;
    end
    xOpt(1) = Time(xstartIndex);
    yOpt(1) = Pressure(xstartIndex);
    
    %% tLeakBeginning
    for i=2:2:length(x)-1
        
       xstartIndex = find(Time < x(i),1,'last');

        if i==length(x)-1 %Top Off Begin
            xOpt(i) = Time(xstartIndex);
            yOpt(i) = Pressure(xstartIndex);
        else
            while Pressure(xstartIndex) <= Pressure(xstartIndex-Frequency) %go left on x axis
                xstartIndex = xstartIndex-Frequency;
            end
            while Pressure(xstartIndex) > Pressure(xstartIndex-3)          %go right on x axis
                xstartIndex = xstartIndex+1;
            end
            xOpt(i) = Time(xstartIndex-2);
            yOpt(i) = Pressure(xstartIndex-2);
        end
    end
    
    %% tLeakEnd
    for i=3:2:length(x)-2
        
     xstartIndex = find(Time < x(i),1,'last');

        while Pressure(xstartIndex) >= Pressure(xstartIndex+Frequency) %go right on x axis
            xstartIndex = xstartIndex+Frequency;
        end
        while Pressure(xstartIndex) < Pressure(xstartIndex+Frequency)  %go left on x axis
            xstartIndex = xstartIndex-1;
        end
        xOpt(i) = Time(xstartIndex+Frequency);
        yOpt(i) = Pressure(xstartIndex+Frequency);
    end
    
    %% tEnd
    xend = length(x);
    xstartIndex = find(Time < x(i),1,'last')

    if Pressure(xstartIndex) <= Pressure(xstartIndex-3)             %go left on x axis
        while Pressure(xstartIndex) <= Pressure(xstartIndex-3)
            xstartIndex = xstartIndex-4*Frequency;
        end
        xOpt(xend) = Time(xstartIndex);
        yOpt(xend) = Pressure(xstartIndex);
    end
    if Pressure(xstartIndex) > Pressure(xstartIndex-2*Frequency)    %go right on x axis
        while Pressure(xstartIndex)<max(Pressure) && xstartIndex<numel(Pressure)
            xstartIndex = xstartIndex+1;
        end
        xOpt(xend) = Time(xstartIndex);
        yOpt(xend) = Pressure(xstartIndex);
    end
    
    
end
% for i=1:length(xOpt)
%     Legendtstart = plot(xOpt(i),yOpt(i),'or');
% end
end