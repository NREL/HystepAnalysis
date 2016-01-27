function [xend, PtargetEnd, xendlastramp, Ptarget1, yendRamp, APRR2] = Plot_Top_Off(table, Ind, MeanAmbTemp, Vdown, Vup, Pzero, APRR, x, y, dPupper, dPlower)

%linear Interpolation for Ptarget1 (transition point for APRR -> APRR Top-Off
Ptarget1up = table(Ind,3); %top-off target, different from non-top off
Ptarget1down = table(Ind+1,3);
Ptarget1 = Ptarget1down + (Ptarget1up - Ptarget1down) * ((MeanAmbTemp - Vdown)/(Vup-Vdown));
%linear Interpolation for Ptarget2
Ptarget2up = table(Ind,4);
Ptarget2down = table(Ind+1,4);
PtargetEnd = Ptarget2down + (Ptarget2up - Ptarget2down) * ((MeanAmbTemp - Vdown)/(Vup-Vdown));
%linear Interpolation for APRR2 (Top-off APRR) 
APRR2up  = table(Ind,5);
APRR2down = table(Ind+1,5);
APRR2Min = APRR2down + (APRR2up - APRR2down) * ((MeanAmbTemp - Vdown)/(Vup-Vdown));
APRR2=APRR2Min/60;

yEndoldRamp=APRR*(x(2)-x(1))+Pzero;

%%% Plot pressure ramp and pressure tolerences
for i=1:2:(length(x)-1)
    if i==1
        xline=[x(i) x(i+1)];
        yline=[y(i) yEndoldRamp];
        yend1=y(i)+dPlower;
        xend1=x(i)+(yend1-Pzero)/APRR;
        xlinestart=[x(i) xend1];
        ylinestart=[Pzero Pzero];
        xlineup=[x(1) x(1)];
        ylineup=[Pzero Pzero+dPupper];
        xendline=[xend1 x(i+1)];
        yendline=[y(i) yEndoldRamp-dPlower];
        yupperTol = yline + dPupper;
        %LegendAPRR = line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
        %LegendAPRRtol = line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
        %line(xendline,yendline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        %line(xlinestart,ylinestart,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        %line(xlineup,ylineup,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        Pc = yline;
        Xc = xline;
        Pupper = [ylineup yupperTol];
        Xupper = [xlineup xline];
        Plower = [ylinestart yendline];
        Xlower = [xlinestart xendline];
    else
        if i==length(x)-1 %for last set
            %%% Plot end of first APRR and top-off
            yendRamp=Ptarget1; %Transistion Pressure
            xendlastramp=x(i)+(yendRamp-yEndoldRamp)/APRR;
            xline=[x(i) xendlastramp];
            yline=[yEndoldRamp yendRamp];
            yupperTol= yline + dPupper;
            ylowerTol = yline - dPlower;
            %line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
            %line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            %line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            Pc = [Pc yline];
            Xc = [Xc xline];
            Pupper = [Pupper yupperTol];
            Xupper = [Xupper xline];
            Plower = [Plower ylowerTol];
            Xlower = [Xlower xline];
            yEndoldRamp=yendRamp;
            
            yendRamp=PtargetEnd;
            xend=xendlastramp+(yendRamp-yEndoldRamp)/APRR2;
            xline=[xendlastramp xend];
            yline=[yEndoldRamp yendRamp];
            yend1=yendRamp-dPupper;
            xend1=xend-(yendRamp-yend1)/APRR2;
            yendline=[yendRamp yendRamp];
            xendline=[xend1 xend];
            xupperTolline=[xendlastramp xend1];
            yupperTolline=[yEndoldRamp+dPupper yendRamp];
            xlineup=[xend xend];
            ylineup=[yendRamp-dPlower yendRamp];
            ylowerTol = yline - dPlower;
            %line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
            %line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            %line(xendline,yendline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            %line(xupperTolline,yupperTolline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            %line(xlineup,ylineup,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            Pc = [Pc yline];
            Xc = [Xc xline];
            Pupper = [Pupper yupperTolline yendline];
            Xupper = [Xupper xupperTolline xendline];
            Plower = [Plower ylowerTol ylineup];
            Xlower = [Xlower xline xlineup];
        else
            yendRamp=APRR*(x(i+1)-x(i))+yEndoldRamp;
            xline=[x(i) x(i+1)];
            yline=[yEndoldRamp yendRamp];
            yupperTol = yline + dPupper;
            ylowerTol = yline - dPlower;
            %line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
            %line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            %line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
            yEndoldRamp=yendRamp;
            Pc = [Pc yline];
            Xc = [Xc xline];
            Pupper = [Pupper yupperTol];
            Xupper = [Xupper xline];
            Plower = [Plower ylowerTol];
            Xlower = [Xlower xline];
        end
    end
end

%Plot Tolerances for Leaks
yEndoldLeak=APRR*(x(2)-x(1))+y(1);
% for i=2:2:(length(x)-1)
%     if i==2
%         yendLeak=APRR*(x(i)-x(i-1))+y(i-1);
%         xline=[x(i) x(i+1)];
%         yline=[yendLeak yendLeak];
%         yupperTol = yline + dPupper;
%         ylowerTol = yline - dPlower;
%         line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
%         line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%         line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%     else
%         if i==length(x)-2
%             yendRamp=PtargetEnd;
%             xend=xendlastramp+(yendRamp-yEndoldRamp)/APRR2;
%             xline=[xendlastramp xend];
%             yline=[yEndoldRamp yendRamp];
%             yend1=yendRamp-dPupper;
%             xend1=xend-(yendRamp-yend1)/APRR2;
%             yendline=[yendRamp yendRamp];
%             xendline=[xend1 xend];
%             xupperTolline=[xendlastramp xend1];
%             yupperTolline=[yEndoldRamp+dPupper yendRamp];
%             xlineup=[xend xend];
%             ylineup=[yendRamp-dPlower yendRamp];
%             ylowerTol = yline - dPlower;
%             line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
%             line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xendline,yendline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xupperTolline,yupperTolline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xlineup,ylineup,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             
%             yendLeak=APRR*(x(i)-x(i-1))+yEndoldLeak;
%             xline=[x(i) x(i+1)];
%             yline=[yendLeak yendLeak];
%             yupperTol = yline + dPupper;
%             ylowerTol = yline - dPlower;
%             line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
%             line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             yEndoldLeak=yendLeak;
%         
%         else
%             yendLeak=APRR*(x(i)-x(i-1))+yEndoldLeak;
%             xline=[x(i) x(i+1)];
%             yline=[yendLeak yendLeak];
%             yupperTol = yline + dPupper;
%             ylowerTol = yline - dPlower;
%             line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2)
%             line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             yEndoldLeak=yendLeak;
%         end
%     end
% end
line(Xc,Pc,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
line(Xupper,Pupper,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
line(Xlower,Plower,'Color',[1 0 0],'LineStyle',':','LineWidth',2);

if length(x) <= 4
    yendRamp=PtargetEnd;
    xend=x(i)+(yendRamp-y(i))/APRR;
    xendlastramp =xend;
    if yline(1,2)<PtargetEnd
        fprintf(2, '\n!!! Leak-Checks or Top Off Fueling missing!!!\n\n')
    end
    
end
end