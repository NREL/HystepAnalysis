%% Written(original) by T. Luetje on 12/22/2014
%Browse for File
%Convert File
%Get right Table
%Import Table
%Extract signals of converted File
%Plot Fill
%Calculate Results
%Show Results

function Plot_Fill_Standard_OneFill_hystep(handles, filename, pathname,sheetnum)
dbstop if error

%% %%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
programPath = '';


%%%% HySTEP Device


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import and Prepare Fill Data File
sourcefile = [pathname filename];

[~,sheets] = xlsfinfo(sourcefile);

[fillData, fillDataHeaders, ~] = xlsread(sourcefile, sheets{sheetnum}); % output is [num, txt, raw]
[overviewData, overviewDataTxt, overviewDataRaw] = xlsread(sourcefile, sheets{1});

%%%% Data Sheet

%%%%% Future improvement - make dynamic and accomodate multiple frequencies
[K,L] = find(strcmp('Channels',overviewDataTxt));
channels = overviewDataRaw{K+sheetnum-1,L};
[I,J] = find(strcmp('wf_increment',overviewDataTxt));
inc = [overviewDataRaw{I(sheetnum-1)+1:I(sheetnum-1)+channels,J(sheetnum-1)}];
inc = unique(inc);
[M,N] = find(strcmp('Length',overviewDataTxt));
varlen = [overviewDataRaw{M(sheetnum-1)+1:M(sheetnum-1)+channels,N(sheetnum-1)}];
varlen = unique(varlen);
inc = sort(inc,'ascend');
varlen = sort(varlen,'descend');
[O,P] = find(strcmp('Group',overviewDataTxt));
testname = overviewDataRaw{O+sheetnum-1,P};
[Q,R] = find(strcmp(testname,overviewDataTxt));
starttime = datenum(datestr(overviewDataTxt{Q(2)+2,13}));

timeStep = min(inc);% [s]
dataFreq = 1/timeStep; % [hz]
t = 0:timeStep:max(varlen)*timeStep-timeStep; %create time vector [s]
%%%%% MULTIPLE TIME VECTORS, MAIN FOR P,T is t
for i = 2:numel(inc)
    eval(['t',num2str(i),' = 0:inc(i):varlen(i)*inc(i)-inc(i);']);
end

for i = 1:length(fillDataHeaders)
    eval([fillDataHeaders{i} ' = fillData(:,strcmp(''' fillDataHeaders{i} ''',fillDataHeaders));']);
end

%clear fillData fillDataHeaders fillDataraw sheets i

Pressure = P_recep;
SFDT = T_recep; %station fuel delivery temperature (closest value)
Time = t;
T_tank1_ave = (T_tank1a+T_tank1b)/2;
T_tank2_ave = (T_tank2a+T_tank2b)/2;
T_tank3_ave = (T_tank3a+T_tank3b)/2;

% Determine which tank is being filled
% Tank pressure is used for determining t0, tleak, and tend.  agnostic to
% number of vessels filled
if TNK1_AV(1) && ~TNK2_AV(1) && ~TNK3_AV(1)
    Ptank = P_tank1;
    Ttank = T_tank1_ave + 273.15;
elseif TNK2_AV(1) && ~TNK1_AV(1) && ~TNK3_AV(1)
    Ptank = P_tank2;
    Ttank = T_tank2_ave + 273.15;
elseif TNK3_AV(1) && ~TNK2_AV(1) && ~TNK1_AV(1)
    Ptank = P_tank3;
    Ttank = T_tank3_ave + 273.15;
elseif TNK1_AV(1) && TNK2_AV(1) && ~TNK3_AV(1)
    Ptank = (P_tank1 + P_tank2)/2;
    Ttank = (T_tank1_ave + T_tank2_ave)/2 + 273.15;
elseif TNK3_AV(1) && TNK2_AV(1) && ~TNK1_AV(1)
    Ptank = (P_tank3 + P_tank2)/2;
    Ttank = (T_tank3_ave + T_tank2_ave)/2 + 273.15;
elseif TNK1_AV(1) && TNK3_AV(1) && ~TNK2_AV(1)
    Ptank = (P_tank1 + P_tank3)/2;
    Ttank = (T_tank1_ave + T_tank3_ave)/2 + 273.15;
else
    Ptank = (P_tank1 + P_tank2 + P_tank3)/3;
    Ttank = (T_tank1_ave + T_tank2_ave + T_tank3_ave)/3 + 273.15;
end
    

%set NWP from GUI input
if get(handles.radiobuttonH70a, 'Value') || get(handles.radiobuttonH70b, 'Value')|| get(handles.radiobuttonH70c, 'Value')
    nwp = 70; %[MPa]
elseif get(handles.radiobuttonH35a, 'Value') || get(handles.radiobuttonH35b, 'Value')
    nwp = 35; %[MPa]
end



%%%% Calculate SOC for confirmation between station calcuation (assuming no
%%%% dispenser offset and on board SOC

%Define Constants
Mmol = 0.00201588;
R = 8.314472;

%Reference Density
Znwp = compressibility(15, nwp);
rhonwp = (nwp*1000000 * Mmol) / (R * (15+273.15) * Znwp); %P [Pa] M / zRT [K]

% if get(handles.radiobuttonComm, 'Value') == 1
%     for i=1:length(MP)
%         if ~isnan(MT(i))
%             Ztemp = compressibility(MT(i)-273.15,MP(i));
%             rhoTemp = (MP(i)*1000000 * Mmol) / (R * (MT(i)) * Ztemp); %Ttank by definition is [K]
%             SOCcalc(i) = rhoTemp/rhonwp*100;% [%]
%             mass_tot(i) = rhoTemp*TV(i);
%         end
%     end
%     mass_delta = diff(decimate(mass_tot,2));
%     mass_delta = [0 mass_delta];
%     flow_rate = mass_delta./(inc(2)*2);
% else
%     valve_tank1 = resample(TNK1_AV(1:length(t2)),length(t),length(t2));
%     valve_tank2 = resample(TNK2_AV(1:length(t2)),length(t),length(t2));
%     valve_tank3 = resample(TNK3_AV(1:length(t2)),length(t),length(t2));
for i=1:length(Ptank)
    if ~isnan(Ttank(i))
        Ztemp = compressibility(Ttank(i)-273.15,Ptank(i));
        rhoTemp = (Ptank(i)*1000000 * Mmol) / (R * (Ttank(i)) * Ztemp); %Ttank by definition is [K]
        SOCcalc(i) = rhoTemp/rhonwp*100;% [%]
        %             mass_tot(i) = rhoTemp*(valve_tank1(i)+valve_tank2(1)+valve_tank3(1))*76;
        mass_tot(i) = rhoTemp*(TNK1_AV(1)+TNK2_AV(1)+TNK3_AV(1))*76;
    end
end
mass_delta = diff(decimate(mass_tot,10));
mass_delta = [0 mass_delta];
flow_rate = mass_delta./(inc(1)*10);
% end

clear Ztemp rhoTemp Znwp rhonwp R Mmol  i
% clear dP1 dP2 dP3

%%  Determine and import specific lookup table to use from GUI Values %%%%
if get(handles.radiobuttonH35a, 'Value') == 1
    if get(handles.radiobuttonComm, 'Value') == 1
        sheet = '35MPa, 2.4-4.2kg (comm)';
        if get(handles.radiobuttonT40, 'Value') == 1
            range = 'C6:M21';
        elseif get(handles.radiobuttonT30, 'Value') == 1
            range = 'C28:M43';
        elseif get(handles.radiobuttonT20, 'Value') == 1
            range = 'C50:M65';
        end
    end
    if get(handles.radiobuttonNonComm, 'Value') == 1
        sheet = '35MPa, 2.4-4.2kg (non-comm)';
        if get(handles.radiobuttonT40, 'Value') == 1
            range = 'C6:M21';
        elseif get(handles.radiobuttonT30, 'Value') == 1
            range = 'C28:M43';
        elseif get(handles.radiobuttonT20, 'Value') == 1
            range = 'C50:M65';
        end
    end
end
if get(handles.radiobuttonH35b, 'Value') == 1
    if get(handles.radiobuttonComm, 'Value') == 1
        sheet = '35MPa, 4.2-6.0kg (comm)';
        if get(handles.radiobuttonT40, 'Value') == 1
            range = 'C6:M21';
            
        elseif get(handles.radiobuttonT30, 'Value') == 1
            range = 'C28:M43';
            
        elseif get(handles.radiobuttonT20, 'Value') == 1
            range = 'C50:M65';
        end
    end
    if get(handles.radiobuttonNonComm, 'Value') == 1
        sheet = '35MPa, 4.2-6.0kg (non-comm)';
        if get(handles.radiobuttonT40, 'Value') == 1
            range = 'C6:M21';
        elseif get(handles.radiobuttonT30, 'Value') == 1
            range = 'C28:M43';
        elseif get(handles.radiobuttonT20, 'Value') == 1
            range = 'C50:M65';
        end
    end
end
if get(handles.radiobuttonH70a, 'Value') == 1
    if get(handles.checkboxCD0C, 'Value') == 1
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 2-4kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C75:S90';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C97:S112';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 2-4kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C75:P90';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C97:P112';
            end
        end
    end
    if get(handles.checkboxCD10C, 'Value') == 1
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 2-4kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C122:S137';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C144:S159';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 2-4kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C122:P137';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C144:P159';
            end
        end
    end
    if  get(handles.checkboxCD0C, 'Value') == 0 && get(handles.checkboxCD10C, 'Value') == 0
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 2-4kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C6:S21';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C28:S43';
            elseif get(handles.radiobuttonT20, 'Value') == 1
                range = 'C50:S65';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 2-4kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C6:P21';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C28:P43';
            elseif get(handles.radiobuttonT20, 'Value') == 1
                range = 'C50:P65';
            end
        end
    end
end
if get(handles.radiobuttonH70b, 'Value') == 1
    if get(handles.checkboxCD0C, 'Value') == 1
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 4-7kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C75:S90';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C97:S112';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 4-7kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C75:P90';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C97:P112';
            end
        end
    end
    if get(handles.checkboxCD10C, 'Value') == 1
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 4-7kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C122:S137';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C144:S159';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 4-7kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C122:P137';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C144:P159';
            end
        end
    end
    if  get(handles.checkboxCD0C, 'Value') == 0 && get(handles.checkboxCD10C, 'Value') == 0
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 4-7kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C6:S21';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C28:S43';
            elseif get(handles.radiobuttonT20, 'Value') == 1
                range = 'C50:S65';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 4-7kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C6:P21';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C28:P43';
            elseif get(handles.radiobuttonT20, 'Value') == 1
                range = 'C50:P65';
            end
        end
    end
end
if get(handles.radiobuttonH70c, 'Value') == 1
    if get(handles.checkboxCD0C, 'Value') == 1
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 7-10kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C75:S90';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C97:S112';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 7-10kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C75:P90';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C97:P112';
            end
        end
    end
    if get(handles.checkboxCD10C, 'Value') == 1
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 7-10kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C122:S137';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C144:S159';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 7-10kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C122:P137';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C144:P159';
            end
        end
    end
    if get(handles.checkboxCD0C, 'Value') == 0 && get(handles.checkboxCD10C, 'Value') == 0
        if get(handles.radiobuttonComm, 'Value') == 1
            sheet = '70MPa, 7-10kg (comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C6:S21';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C28:S43';
            elseif get(handles.radiobuttonT20, 'Value') == 1
                range = 'C50:S65';
            end
        end
        if get(handles.radiobuttonNonComm, 'Value') == 1
            sheet = '70MPa, 7-10kg (non-comm)';
            if get(handles.radiobuttonT40, 'Value') == 1
                range = 'C6:P21';
            elseif get(handles.radiobuttonT30, 'Value') == 1
                range = 'C28:P43';
            elseif get(handles.radiobuttonT20, 'Value') == 1
                range = 'C50:P65';
            end
        end
    end
end

LTfilepath = ([ programPath '\Look up Tables Standard.xlsx']);
[LTdata, LTtxt, LTraw] = xlsread(LTfilepath, sheet, range);


%% get graphical user input to determine tstart, tleak, and tend

%Create Figure with Fill
fig = figure('Color', 'white','visible','on');

LegendFillGraph = plot(Time,Ptank);

% set(fig,'units','normalized','outerposition',[0 0 1 1])
xlabel('Time [sec]','Fontsize',12);
ylabel('Pressure [MPa]','Fontsize',12);
title('Please Mark Start of Fill, Start & End of Leaks & End of Fill and Press Return','Color','r','Fontsize',16);
hold on

%crop fill for user input
% PstartIndex = find(Ptank == min(Ptank),1,'last');
% PEndIndex = find(Ptank == max(Ptank),1,'first');
% axis([(PstartIndex/dataFreq)-25 PEndIndex/dataFreq+25 0 85])

%Get input for tStart, tLeak, tEnd
[x,y] = ginput;

% Function for finding exact Values (optimizing user input)
[x, y]=Exact_Values_FRR(x, y, Ptank, Time, dataFreq);
%[x, y, Legendtstart]=Exact_Values_FRR(x, y, Ptank, Time, dataFreq);

idxStart = find(Time == x(1)); %150830 only valid for main time freq
idxEnd = find(Time == x(end));

T_Amb = T_amb(idxStart - dataFreq); %take 1s before start of fill (assuming data starts near when dispenser calculates table values)


%% Get needed Signals & Values
Pzero = y(1);
dPupper = 7.0;
dPlower = 2.5;

% find Ambient Temperature in Table
tableRow = find(T_Amb < LTdata(:,1),1,'last');
tableColumn = find(Pzero > LTdata(1,:),1,'last');

Vup = LTdata(tableRow,1); %upper ambient temperature boundary
Vdown = LTdata(tableRow+1,1);%lower ambient temperature boundary
APRRup = LTdata(tableRow,2);
APRRdown = LTdata(tableRow+1,2);

%linear Interpolation for APPR
APRRMin = APRRdown + (APRRup - APRRdown) * ((T_Amb - Vdown) / (Vup - Vdown));
APRR = APRRMin/60; %APPR in MPa/s.  Used for calculating tolerences



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Target Pressure with Interpolation depending on Pzero %%%%%%%%%%%%%%%%

% ColumnIndex = max(logIndex);
Pdown = LTdata(1,tableColumn);       % Pressure Support Value which is fewer than the actual Ambient Temperature Tamb
Pup = LTdata(1,tableColumn+1);  

% Pressure Support Value which is greater than the actual Ambient Temperature Tamb
luValue = LTdata(tableRow,tableColumn);         % left upper value in lookup table
llValue = LTdata(tableRow+1,tableColumn);       % left lower value in lookup table
ruValue = LTdata(tableRow,tableColumn+1);       % right upper value in lookup table
rlValue = LTdata(tableRow+1,tableColumn+1);     % right lower value in lookup table

%%%%%%%%%%%%%%%%%%%
% Cases: "70MPa, -40 to -33°C, 2-4kg (SOC = 100%)" and "70MPa, -26 to % -17,5°C, 7-10kg" can not yet be calculated properly in Case Tamb > 45 degree



if (isnan(luValue) == 1 || isnan(llValue) == 1 || isnan(ruValue) == 1 || isnan(rlValue) == 1) % Checks if any of the four values are NaN
    
    % For ambient temperatures, where at least one interpolation point for Fueling Target Pressure states “no top-off”
    
    if (strcmp(LTraw(tableRow+2,tableColumn+1),'see Top-Off') && ~ strcmp(LTraw(tableRow+3,tableColumn+1),'see Top-Off')) % Checks if left upper is, and if left lower vaule is not "see Top Off", no top off shall be applied
        
        luValue = LTdata(RowIndex,3);          % replaces left upper value in Top Off Value as discussed in SAE J2601 page 92
        ruValue = LTdata(RowIndex,3);          % replaces right upper value in Top Off Value as discussed in SAE J2601 page 92
        
        Ptargetdown = llValue + (luValue - llValue) * ((T_Amb - Vdown)/(Vup - Vdown));           % linear Interpolation between left upper and left lower Value
        Ptargetup = rlValue + (ruValue - rlValue) * ((T_Amb - Vdown)/(Vup - Vdown));             % linear Interpolation between right upper and right lower Value
        PtargetEnd = Ptargetdown + (Ptargetup - Ptargetdown) * ((Pup - Pzero)/(Pup - Pdown));   % linear Interpolation between the two interim results above
        
        flagTopOff = false;
        
    elseif (strcmp(LTraw(tableRow+2,tableColumn+1),'see Top-Off') && strcmp(LTraw(tableRow+3,tableColumn+1),'see Top-Off'))   % Top Off Fueling is the Case
        
        flagTopOff = true;
        
    end
    
else
    
    % linear Interpolation (between 4 Values out of Table) for PTDown, PTup and then PT
    
    Ptargetdown = llValue + (luValue - llValue) * ((T_Amb - Vdown)/(Vup - Vdown));           % linear Interpolation between left upper and left lower Value
    Ptargetup = rlValue + (ruValue - rlValue) * ((T_Amb - Vdown)/(Vup - Vdown));             % linear Interpolation between right upper and right lower Value
    PtargetEnd = Ptargetdown + (Ptargetup - Ptargetdown) * ((Pup - Pzero)/(Pup - Pdown));   % linear Interpolation between the two interim results above
    
    flagTopOff = false;
    
end

% %% Plot Tolerances for Top Off Fueling
if flagTopOff
    
%     [xend, PtargetEnd, xendlastramp, Ptarget1, yendRamp, LegendAPRR, LegendAPRRtol, APRR2] = ...               %CALL FUNCTION (Plot Top Off)
%         Plot_Top_Off(LTdata, tableRow, T_Amb, tAmbDown, tAmbUp, Pzero, APRR, x, y, dPupper, dPlower) ;
    [xend, PtargetEnd, xendlastramp, Ptarget1, yendRamp, APRR2] = ...               %CALL FUNCTION (Plot Top Off)
        Plot_Top_Off(LTdata, tableRow, T_Amb, Vdown, Vup, Pzero, APRR, x, y, dPupper, dPlower) ;
    
else
    
    %% PTarget Interpolations for Ptarget limit Values to TOF (See Tables)
    
    %Get P_Target from Input
    Index = find(Pzero>LTdata(1,:),1,'last');
    Pdown = LTdata(1,Index);
    Pup = LTdata(1,Index+1);
    
    % Ptarget Top Off 70MPa 2-4kg (comm) T40
    if  strcmp(sheet,'70MPa, 2-4kg (comm)')==1 && strcmp(range,'C6:S21')==1 && (T_Amb<10 && T_Amb>0 || T_Amb>45 && T_Amb<50) && Pzero>=0.5 && Pzero<=5;
        if T_Amb > 45
            PtargetEnd=LTdata(tableRow+1,Index+1);
        else
            Ptargetup=LTdata(tableRow,Index+1);
            Ptargetdown=LTdata(tableRow,Index);
            PtargetEnd = Ptargetdown + (Ptargetup - Ptargetdown) * ((Pdown - Pzero) / (Pup-Pdown));
        end
    else
        
        if  (strcmp(sheet,'70MPa, 2-4kg (comm)')==1 && strcmp(range,'C28:S43')==1 && T_Amb<-30 && Pzero>=0.5 && Pzero<=5) || ... % Ptarget Top Off 70MPa 2-4kg (comm) T30
                (strcmp(sheet,'70MPa, 2-4kg (comm)')==1 && (strcmp(range,'C75:S90')==1 || strcmp(range,'C122:S137')==1) && T_Amb<10 && T_Amb>0 && Pzero>=0.5 && Pzero<=5) || ... % Ptarget Top Off 70MPa 2-4kg (comm) T40 CD0
                (strcmp(sheet,'70MPa, 2-4kg (comm)')==1 && (strcmp(range,'C97:S112')==1 || strcmp(range,'C144:S159')==1) && T_Amb<-30 && T_Amb>-40 && Pzero>=0.5 && Pzero<=5) || ... % Ptarget Top Off 70MPa 2-4kg (comm) T30 CD0
                (strcmp(sheet,'70MPa, 4-7kg (comm)')==1 && (strcmp(range,'C6:S21')==1 || strcmp(range,'C75:S90')==1 || strcmp(range,'C122:S137')==1) && T_Amb<10 && T_Amb>0 && Pzero>=0.5 && Pzero<=5) || ... % Ptarget Top Off 70MPa 4-7kg (comm) T40
                (strcmp(sheet,'70MPa, 7-10kg (comm)')==1 && strcmp(range,'C6:S21')==1 && T_Amb<20 && T_Amb>10 && Pzero>=0.5 && Pzero<=5) || ... %Ptarget Top Off 70MPa 7-10kg (comm) T40
                (strcmp(sheet,'70MPa, 7-10kg (comm)')==1 && strcmp(range,'C50:S65')==1 && T_Amb>45 && Pzero>=0.5 && Pzero<=5) || ... % Top Off 70MPa 7-10kg (comm) T20
                (strcmp(sheet,'70MPa, 7-10kg (comm)')==1 && strcmp(range,'C75:S90')==1 && T_Amb<25 && T_Amb>20 && Pzero>=0.5 && Pzero<=5) || ... % Ptarget Top Off 70MPa 7-10kg (comm) T40 & CD0
                (strcmp(sheet,'70MPa, 7-10kg (comm)')==1 && strcmp(range,'C122:S137')==1 && T_Amb<30 && T_Amb>25 && Pzero>=0.5 && Pzero<=5) % Ptarget Top Off 70MPa 7-10kg (comm) T40 & CD-10
            Ptargetup=LTdata(tableRow,Index+1);
            Ptargetdown=LTdata(tableRow,Index);
            PtargetEnd = Ptargetdown + (Ptargetup - Ptargetdown) * ((Pdown - Pzero)/(Pup - Pdown));
        else
            %Get P_Target from Input
            PtargetAmbdown = LTdata(tableRow,Index);
            PtargetAmbup = LTdata(tableRow+1,Index);
            PtargetPdown = LTdata(tableRow,Index+1);
            PtargetPup = LTdata(tableRow+1,Index+1);
            
            %linear Interpolation (between 4 Values out of Table) for PTDown, PTup and then PT
            Ptargetdown = PtargetAmbup + (PtargetAmbdown - PtargetAmbup) * ((T_Amb - Vdown)/(Vup-Vdown));
            Ptargetup = PtargetPup + (PtargetPdown - PtargetPup) * ((T_Amb - Vdown)/(Vup-Vdown));
            PtargetEnd = Ptargetdown + (Ptargetup - Ptargetdown) * ((y(1)-Pdown)/(Pup-Pdown));
        end
    end
    
    
    
    
    %% Plot Tolerances without Top Off Fueling
    %TIR stuff **********************************************************
    %PtargetEnd = 72.2;
    %APRR = 9.14/60;
    %***********************************************************************
    %Plot Tolerances for Ramps
    if length(x)==2
        yendRamp=PtargetEnd;
        xend=x(1)+(yendRamp-y(1))/APRR;
        xline=[x(1) xend];
        yline=[y(1) yendRamp];
        yend2=yendRamp-dPupper;
        xend2=xend-(yendRamp-yend2)/APRR;
        yend1=y(1)+dPlower;
        xend1=x(1)+(yend1-Pzero)/APRR;
        xlinestart=[x(1) xend1];
        ylinestart=[Pzero Pzero];
        xlineup1=[x(1) x(1)];
        ylineup1=[Pzero Pzero+dPupper];
        xstartline=[xend1 xend];
        ystartline=[y(1) yendRamp-dPlower];
        xendline=[x(1) xend2];
        yendline=[Pzero+dPupper yendRamp];
        xlineup2=[xend xend];
        ylineup2=[yendRamp-dPlower PtargetEnd];
        xlineEnd=[xend2 xend];
        ylineEnd=[yendRamp yendRamp];
        %         LegendAPRR = line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
        %         LegendAPRRtol = line(xstartline,ystartline,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
        %         line(xlinestart,ylinestart,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        %         line(xlineup1,ylineup1,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        %         line(xlineup2,ylineup2,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        %         line(xendline,yendline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        %         line(xlineEnd,ylineEnd,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
        Pc = yline;
        Xc = xline;
        Pupper = [ylineup1 yendline ylineEnd];
        Xupper = [xlineup1 xendline xlineEnd];
        Plower = [ylinestart ystartline ylineup2];
        Xlower = [xlinestart xstartline xlineup2];
    else
        yEndoldRamp=APRR*(x(2)-x(1))+Pzero;
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
                %             LegendAPRR = line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
                %             LegendAPRRtol = line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
                %             line(xendline,yendline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                %             line(xlinestart,ylinestart,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                %             line(xlineup,ylineup,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                Pc = yline;
                Xc = xline;
                Pupper = [ylineup yupperTol];
                Xupper = [xlineup xline];
                Plower = [ylinestart yendline];
                Xlower = [xlinestart xendline];
            else
                if i==length(x)-1
                    yendRamp=PtargetEnd;
                    xend=x(i)+(yendRamp-yEndoldRamp)/APRR;
                    xline=[x(i) xend];
                    yline=[yEndoldRamp yendRamp];
                    yend1=yendRamp-dPupper;
                    xend1=xend-(yendRamp-yend1)/APRR;
                    yendline=[yendRamp yendRamp];
                    xendline=[xend1 xend];
                    xupperTolline=[x(i) xend1];
                    yupperTolline=[yEndoldRamp+dPupper yendRamp];
                    xlineup=[xend xend];
                    ylineup=[yendRamp-dPlower yendRamp];
                    ylowerTol = yline - dPlower;
                    %                 LegendAPRR =     line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
                    %                 line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                    %                 line(xendline,yendline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                    %                 line(xupperTolline,yupperTolline,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                    %                 line(xlineup,ylineup,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
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
                    %                 LegendAPRR = line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
                    %                 line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
                    %                 line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
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
    end
    
    %     LegendAPRR = line(Xc,Pc,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
    %     LegendAPRRtol = line(Xupper,Pupper,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
    line(Xc,Pc,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
    % This is for TIR stuff
%*************************************************************************************************    
%     Pcupper(1) = Pc(1);
%     Pclower(1) = Pc(1);
%     for i = 2:length(Pc)
%         Pcupper(i) = Pcupper(i-1)+1.1*(Pc(i)-Pc(i-1));
%         Pclower(i) = Pclower(i-1)+0.9*(Pc(i)-Pc(i-1));
%     end
%     line(Xc,Pcupper,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
%     line(Xc,Pclower,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
%*************************************************************************************************       
    line(Xupper,Pupper,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
    line(Xlower,Plower,'Color',[1 0 0],'LineStyle',':','LineWidth',2);
    
    %Plot Tolerances for Leaks
    yEndoldLeak=APRR*(x(2)-x(1))+y(1);
%     for i=2:2:(length(x)-1)
%         if i==2
%             yendLeak=APRR*(x(i)-x(i-1))+y(i-1);
%             xline=[x(i) x(i+1)];
%             yline=[yendLeak yendLeak];
%             yupperTol = yline + dPupper;
%             ylowerTol = yline - dPlower;
%             LegendAPRR = line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
%             line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             d=i;        % for verification in Line 715 (if xend == x(d+1);)
%         else
%             yendLeak=APRR*(x(i)-x(i-1))+yEndoldLeak;
%             xline=[x(i) x(i+1)];
%             yline=[yendLeak yendLeak];
%             yupperTol = yline + dPupper;
%             ylowerTol = yline - dPlower;
%             LegendAPRR = line(xline,yline,'Color',[0 1 0],'LineStyle',':','LineWidth',2);
%             line(xline,yupperTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             line(xline,ylowerTol,'Color',[1 0 0],'LineStyle',':','LineWidth',2)
%             yEndoldLeak=yendLeak;
%             d=i;        % for verification in Line 715 (if xend == x(d+1);)
%         end
%     end
end

%% Format & Save Plot

if exist ('xend','var');
    if x(end)<=xend
        xlimit=xend+10;
    else xlimit=x(end)+60;
    end
else xlimit=x(end)+60;
end
if exist('yendRamp','var')
    yend=yendRamp;
end
if exist ('yend','var');
    if y(end)<=yend
        ylimit=yend+5;
    else ylimit=y(end)+5;
    end
else ylimit=y(end)+5;
end


grid on
if get(handles.checkboxCD0C,'Value')==1;
    ColdDisp= '0°C';
else
    if get(handles.checkboxCD10C,'Value')==1;
        ColdDisp= '-10°C';
    else
        ColdDisp= 'No';
    end
end

% AmbLinex=[x(1) x(end)];
% AmbLiney=[T_Amb T_Amb];
% LegendAmbient = line(AmbLinex,AmbLiney,'Color',[0 0 0]);

%% Plot SFDT and Boundaries
upBoundx=[x(1) x(end)];
downBoundx=[x(1) x(end)];
ylimitdown=-45;

if get(handles.radiobuttonT40,'Value')==1;
    DispTemp= 'T40 (-40°C to -33°C)';
    upBoundy=[-33 -33];
    downBoundy=[-40 -40];
    
elseif get(handles.radiobuttonT30,'Value')==1;
    DispTemp= 'T30 (-33°C to -26°C)';
    upBoundy=[-26 -26];
    downBoundy=[-33 -33];
    
elseif get(handles.radiobuttonT20,'Value')==1;
    DispTemp= 'T20 (-26°C to -17.5°C)';
    upBoundy=[-17.5 -17.5];
    downBoundy=[-26 -26];
    %downBoundy=[-22.5 -22.5]; %TIR stuff
end
delete(LegendFillGraph);
%delete(Legendtstart);
% disp_data = input(['Have NREL dispenser data? y/n [y]:'],'s');
% if isempty(disp_data)
%     disp_data = 'y';
% end
% if disp_data == 'y'
%     text_file_convert
%     rampstart = starttime+x(1)/(60*60*24);
%     figure;plot((D(1).HosePressureH70Time),D(1).HosePressureH70/145.04,'LineWidth',2);
%     axis([starttime+0.95/24 starttime+1.05/24 0 70]);
%     title('Please Zoom in, Mark Start of Fill and Press Return','Color','r','Fontsize',16);
%     [xx yy] = ginput;
%     close
%     kk = find(D.HoseTemperatureH70Time > xx & D.HoseTemperatureH70Time < (xx + t(end)/(60*60*24)));
%     %    LegendDpressure = plot((D(1).HosePressureH70Time(kk) - D(1).HosePressureH70Time(kk(1)))*60*60*24,D(1).HosePressureH70(kk)/145.04,'LineWidth',2);
%     %    LegendDtemp = plot((D(1).HoseTemperatureH70Time(kk) - D(1).HoseTemperatureH70Time(kk(1)))*60*60*24,D(1).HoseTemperatureH70(kk),'LineWidth',2);
%     plot((D(1).HosePressureH70Time(kk) - D(1).HosePressureH70Time(kk(1)))*60*60*24 + x(1),D(1).HosePressureH70(kk)/145.04,'LineWidth',2);
%     %plot((D(1).HosePressureH70Time-(rampstart+0.041833514347672))*60*60*24+x(1),D(1).HosePressureH70/145.04,'LineWidth',2);
% end
% LegendHpressure = plot(Time(idxStart:idxEnd), P_recep(idxStart:idxEnd),'c','LineWidth',2);
% LegendT1pressure = plot(Time(idxStart:idxEnd), P_tank1(idxStart:idxEnd),'LineWidth',2);
% LegendT2pressure = plot(Time(idxStart:idxEnd), P_tank2(idxStart:idxEnd),'LineWidth',2);
% LegendT3pressure = plot(Time(idxStart:idxEnd), P_tank3(idxStart:idxEnd),'LineWidth',2);
% 
% plot(Time(idxStart:idxEnd), P_recep(idxStart:idxEnd),'c','LineWidth',2);
% plot(Time(idxStart:idxEnd), P_tank1(idxStart:idxEnd),'LineWidth',2);
% plot(Time(idxStart:idxEnd), P_tank2(idxStart:idxEnd),'LineWidth',2);
% plot(Time(idxStart:idxEnd), P_tank3(idxStart:idxEnd),'LineWidth',2);
 plot(Time, P_recep,'c','LineWidth',2);
 plot(Time, P_tank1,'LineWidth',2);
 plot(Time, P_tank2,'LineWidth',2);
 plot(Time, P_tank3,'LineWidth',2);
 if nwp == 35
     SOC_1 = SOC_1*40.2/24;
     SOC_2 = SOC_2*40.2/24;
     SOC_3 = SOC_3*40.2/24;
 end
if exist('t3','var')
    %     LegendT1SOC = plot(t3, SOC_1(1:length(t3)),'LineWidth',2);
    %     LegendT2SOC = plot(t3, SOC_2(1:length(t3)),'LineWidth',2);
    %     LegendT3SOC = plot(t3, SOC_3(1:length(t3)),'LineWidth',2);
    plot(t3, SOC_1(1:length(t3)),'LineWidth',2);
    plot(t3, SOC_2(1:length(t3)),'LineWidth',2);
    plot(t3, SOC_3(1:length(t3)),'LineWidth',2);
else
    %     LegendT1SOC = plot(t2, SOC_1(1:length(t2)),'LineWidth',2);
    %     LegendT2SOC = plot(t2, SOC_2(1:length(t2)),'LineWidth',2);
    %     LegendT3SOC = plot(t2, SOC_3(1:length(t2)),'LineWidth',2);
    plot(t2, SOC_1(1:length(t2)),'LineWidth',2);
    plot(t2, SOC_2(1:length(t2)),'LineWidth',2);
    plot(t2, SOC_3(1:length(t2)),'LineWidth',2);
end


% Show when HALT or ABORT signals are sent
if exist('FC','var')
    k = find(FC==3);
    if ~isempty(k)
        plot(t2(k(1)),0,'*k')
        text(t2(k(1))+1,0,'ABORT','FontSize',8)
        line([t2(k(1)) t2(k(1))],[-40 80],'Color','k','LineStyle',':','LineWidth',1)
    end
    k = find(FC==2);
    if ~isempty(k)
        plot(t2(k(1)),0,'*k')
        text(t2(k(1))+1,0,'HALT','FontSize',8)
        line([t2(k(1)) t2(k(1))],[-40 80],'Color','k','LineStyle',':','LineWidth',1)
        l = find(FC(k(end):end)==0);
        if ~isempty(l)
            plot(t2(l(1)+k(end)),0,'*k')
            text(t2(l(1)+k(end))+1,0,'RESUME','FontSize',8)
            line([t2(l(1)+k(end)) t2(l(1)+k(end))],[-40 80],'Color','k','LineStyle',':','LineWidth',1)
        end
    end
end

filltime = t(idxEnd) - t(idxStart);

% plot bars to show leak checks
if length(x)>2
    for i = 1:(length(x)-2)/2
        colormap(gray)
        surf([x(2*i) x(2*i+1)],[0 100],[0.5 0.5;0.5 0.5],'FaceAlpha',0.5)
        filltime = filltime - (x(2*i+1)-x(2*i));
    end
end
%LegendSFDTtol = line(upBoundx,upBoundy,'Color',[1 0 0],'LineStyle','--');
mass_prefill = mass_tot(idxStart) - mass_tot(1);
mass_filled = mass_tot(idxEnd) - mass_tot(idxStart);
text(0,90,{['Startup Mass of H2 = ' num2str(round(mass_prefill,3,'significant')) ' g'];['Mass of H2 Filled = ' num2str(round(mass_filled/1000,3,'significant')) ' kg'];['Fill Time = ' num2str(filltime) ' sec']},'EdgeColor',[0 0 0],'BackgroundColor',[1 1 1],'FontSize',12)

Filename = strrep(filename, '_', '\_' );
title({['Test: ' sheets{sheetnum}];['Version: Standard ' LTraw{1,1}(1:7) ' for ' sheet]; ['APRR = ' num2str(round(APRR*60,3,'significant')) ' MPa/min' ', Target Pressure = ' num2str(round(PtargetEnd,3,'significant')) ' MPa']},'Color','k','Interpreter','none');
% legend([Legendtstart, LegendAPRR, LegendAPRRtol, LegendHpressure, LegendT1pressure, LegendT2pressure, LegendT3pressure, LegendSFDT, LegendT1temp, LegendT2temp, LegendT3temp,LegendT1SOC,LegendT2SOC,LegendT3SOC, LegendSFDTtol,LegendAmbient],...
%     't_S_t_a_r_t, t_L_e_a_k and t_E_n_d','APRR_t_a_r_g_e_t','APRR Tol','P_r_e_c','P_t_a_n_k_1','P_t_a_n_k_2','P_t_a_n_k_3','T_r_e_c','T_t_a_n_k_1','T_t_a_n_k_2','T_t_a_n_k_3','SOC_1','SOC_2','SOC_3','SFDT Tol','T_a_m_b','Location','northwest');
% %       'Vehicle Pressure','t_S_t_a_r_t, t_L_e_a_k and t_E_n_d','APRR_t_a_r_g_e_t','APRR Tolerances','Ambient Temperature','Station Hose Pressure','Station Fuel Delivery Temperature', 'Station Fuel Delivery Temperature Tolerances','Location','northwest');
% 

MinIdx={'-1 Min','0 Min' '1 Min' '2 Min' '3 Min' '4 Min' '5 Min' '6 Min' '7 Min' '8 Min' '9 Min' '10 Min' '11 Min' '12 Min' '13 Min' '14 Min' '15 Min' '16 Min' '17 Min' '18 Min' '19 Min' '20 Min'};
axis([x(1)-60 x(end)+60 0 100])          %axis([xmin xmax ymin ymax])
set(gca,'XTick',x(1)-60:60:x(end)+60)         %set(gca, 'XTick', [debut:step:fin])
set(gca,'XTickLabel',MinIdx)
xlabel('Time [min]','Fontsize',12);
ylabel('Pressure [MPa]/ SOC [%]','Fontsize',12);
legend('APRR','Upper Limit','Lower Limit','Receptacle Pressure','Tank 1 Pressure','Tank 2 Pressure','Tank 3 Pressure','Tank 1 SOC','Tank 2 SOC','Tank 3 SOC','Location','EastOutside')
% if disp_data == 'y'
%     legend('APRR','Upper Limit','Lower Limit','Dispenser Pressure','Receptacle Pressure','Tank 1 Pressure','Tank 2 Pressure','Tank 3 Pressure','Tank 1 SOC','Tank 2 SOC','Tank 3 SOC','Location','EastOutside')
% else
%     legend('APRR','Upper Limit','Lower Limit','Receptacle Pressure','Tank 1 Pressure','Tank 2 Pressure','Tank 3 Pressure','Tank 1 SOC','Tank 2 SOC','Tank 3 SOC','Location','EastOutside')
% end
%axis auto
IdxSuffix = strfind(filename,'.'); %find filename suffix
filename(IdxSuffix:end) = [];

saveas(fig,[pathname sheets{sheetnum} ' Pressure'])
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 11 8.5];
fig.PaperPositionMode = 'manual';
% print('5by3DimensionsFigure','-dpng','-r0')
print(fig,[pathname sheets{sheetnum} ' Pressure'],'-dpng','-r0')
% saveas(fig,[pathname filename ' J2601 Performance Validation .png'])
%saveas(fig,[pathname filename ' J2601 Performance Validation ' sheet ]);         %Save Plot
%saveas(fig,[pathname filename ' J2601 Performance Validation ' sheet '.png']);

% New figure for temperatures
fig2 = figure('Color','white');hold on
line(upBoundx,upBoundy,'Color',[1 0 0],'LineStyle','--');
line(downBoundx,downBoundy,'Color',[1 0 0],'LineStyle','--');
% LegendSFDT = plot(Time(idxStart:idxEnd), SFDT(idxStart:idxEnd),'m','LineWidth',2);
% LegendT1temp = plot(Time(idxStart:idxEnd), T_tank1_ave(idxStart:idxEnd),'LineWidth',2);
% LegendT2temp = plot(Time(idxStart:idxEnd), T_tank2_ave(idxStart:idxEnd),'LineWidth',2);
% LegendT3temp = plot(Time(idxStart:idxEnd), T_tank3_ave(idxStart:idxEnd),'LineWidth',2);
% plot(Time(idxStart:idxEnd), SFDT(idxStart:idxEnd),'m','LineWidth',2);
% plot(Time(idxStart:idxEnd), T_tank1_ave(idxStart:idxEnd),'LineWidth',2);
% plot(Time(idxStart:idxEnd), T_tank2_ave(idxStart:idxEnd),'LineWidth',2);
% plot(Time(idxStart:idxEnd), T_tank3_ave(idxStart:idxEnd),'LineWidth',2);
% plot(Time(idxStart:idxEnd), T_amb(idxStart:idxEnd),'LineWidth',2);

plot(Time, SFDT,'m','LineWidth',2);
plot(Time, T_tank1_ave,'LineWidth',2);
plot(Time, T_tank2_ave,'LineWidth',2);
plot(Time, T_tank3_ave,'LineWidth',2);
plot(Time, T_amb,'LineWidth',2);
% if disp_data == 'y'
% plot((D(1).HoseTemperatureH70Time(kk) - D(1).HoseTemperatureH70Time(kk(1)))*60*60*24 + x(1),D(1).HoseTemperatureH70(kk),'LineWidth',2);
% end

% if get(handles.radiobuttonComm, 'Value') == 1
%     plot(decimate(t2,2),flow_rate,'LineWidth',2)
% else
%     plot(t2,flow_rate,'LineWidth',2)
% end
plot(decimate(t,10),flow_rate,'LineWidth',2)
% plot bars to show leak checks
if length(x)>2
    for i = 1:(length(x)-2)/2
        colormap(gray)
        surf([x(2*i) x(2*i+1)],[-45 85],[0.3 0.3;0.3 0.3],'FaceAlpha',0.5)
    end
end
title({['Test: ' sheets{sheetnum}];['Version: Standard ' LTraw{1,1}(1:7) ' for ' sheet];},'Color','k','Interpreter','none','Fontsize',16);
MinIdx={'-1 Min','0 Min' '1 Min' '2 Min' '3 Min' '4 Min' '5 Min' '6 Min' '7 Min' '8 Min' '9 Min' '10 Min' '11 Min' '12 Min' '13 Min' '14 Min' '15 Min' '16 Min' '17 Min' '18 Min' '19 Min' '20 Min'};
axis([x(1)-60 x(end)+60 -45 85])          %axis([xmin xmax ymin ymax])
set(gca,'XTick',x(1)-60:60:x(end)+60)         %set(gca, 'XTick', [debut:step:fin])
set(gca,'XTickLabel',MinIdx)
xlabel('Time [min]','Fontsize',12);
ylabel('Temperature [C]/Mass Flow [g/sec]','Fontsize',12);
grid on
legend('Upper Limit','Lower Limit','Receptacle','Tank 1','Tank 2','Tank 3','Ambient','Mass Flow','Location','EastOutside')
% if disp_data == 'y'
%     legend('Upper Limit','Lower Limit','Receptacle','Tank 1','Tank 2','Tank 3','Ambient','Dispenser','Mass Flow','Location','EastOutside')
% else
%     legend('Upper Limit','Lower Limit','Receptacle','Tank 1','Tank 2','Tank 3','Ambient','Mass Flow','Location','EastOutside')
% end

saveas(fig2,[pathname sheets{sheetnum} ' Temp'])
fig2.PaperUnits = 'inches';
fig2.PaperPosition = [0 0 11 8.5];
fig2.PaperPositionMode = 'manual';
% print('5by3DimensionsFigure','-dpng','-r0')
print(fig2,[pathname sheets{sheetnum} ' Temp'],'-dpng','-r0')

% %% Calculate Results
% Commented the rest of this out for now. Not sure if we want to use it.
% TAJ 09/22/2015
% 
% %   Percentage Pressure Outside Tolerance
% n=zeros(length(x),1);
% for i=1:length(x)
%     IndexArray=x(i)>Time;                                         %Array of Time Values <x(1,1)
%     x1 = numel(IndexArray(IndexArray))-1;
%     n(i)=x1;
% end
% 
% yAr=zeros(max(n),1);                                            %Empty Array for APRR_target-Values
% Pnew=Pzero;
% if exist('APRR2','var')                                         %If Top Off Fueling exists
%     for i=1:2:length(x)-3
%         for k=n(i):n(i+1)
%             y1=APRR*timeStep+Pnew;
%             yAr(k,1)=y1;
%             Pnew=y1;
%         end
%         for k=n(i+1):n(i+2)
%             yAr(k,1)=max(yAr);
%         end
%     end
%     Pnew1=Pnew;
%     if xendlastramp>x(end)
%         for k=n(length(x)-2):n(length(x))
%             y1=APRR*timeStep+Pnew1;
%             yAr(k,1)=y1;
%             Pnew1=y1;
%         end
%     else
%         for k=n(length(x)-2,1):n(length(x)-1,1)                     %still to improve for begin of TOF before x(end,1)
%             y1=APRR*timeStep+Pnew1;
%             yAr(k,1)=y1;
%             Pnew1=y1;
%         end
%         for k=n(length(x)-1,1):n(length(x))
%             y1=APRR2*timeStep+Pnew1;
%             yAr(k,1)=y1;
%             Pnew1=y1;
%         end
%     end
% else
%     for i=1:2:length(x)-3
%         for k=n(i,1):n(i+1,1)-1
%             y1=APRR*timeStep+Pnew;
%             yAr(k,1)=y1;
%             Pnew=y1;
%         end
%         for k=n(i+1,1):n(i+2,1)
%             yAr(k,1)=max(yAr);
%         end
%     end
%     for k=n(length(x)-1,1):n(length(x))
%         y1=APRR*timeStep+Pnew;
%         yAr(k,1)=y1;
%         Pnew=y1;
%     end
% end
% countout=0;
% countin=0;
% for i=n(1,1):n(end,1)
%     if Pressure(i)>yAr(i)+7.5||Pressure(i)<yAr(i)-2.5
%         countout=countout+1;
%     else
%         countin=countin+1;
%     end
% end
% result=countout/(countin+countout)*100;
% if result>0
%     result=round(result);
% end
% PercentageTol = 0;
% if result == PercentageTol
%     PassPercentage=' |Yes |';
% else if result > PercentageTol
%         PassPercentage='|No  |';
%     end
%     if result < PercentageTol
%         fprintf(2, '\n!!! Unreasonable Value for Pressure Tolerance Result!!!\n\n')
%     end
% end
% 
% %Tolerance for Filltime
% FillTime=x(end,1);
% if exist('APRR2','var') && y(end,1)>Ptarget1;
%     %Top off Fill Time Tolerance
%     Pnew=APRR*(x(end,1)-x(end-2,1))+Pnew;
%     xendAPRRtarget=x(end,1)+(y(end,1)-Pnew)/APRR2;
%     yend1=Pnew+dPlower;
%     xend1=xendAPRRtarget+(yend1-Pnew)/APRR2;
%     yend2=y(end,1)-dPupper;
%     xend2=xendAPRRtarget-(y(end,1)-yend2)/APRR2;
%     if FillTime>=xend2 && FillTime<=xend1
%         PassFillTime='Yes |';
%     else
%         PassFillTime='No  |';
%     end
%     xendAPRRtarget=num2str(xendAPRRtarget);
% else
%     %Ptarget1 FillTime Tolerance
%     if exist('APRR2','var') && xendlastramp>x(end,1)
%         Pnew=APRR*(x(end,1)-x(end-2,1))+Pnew;
%     else
%         if exist('APRR2','var') && xendlastramp<=x(end,1)
%             Pnew=APRR*(xendlastramp-x(end-2,1))+Pnew;
%             Pnew=APRR2*(x(end,1)-xendlastramp)+Pnew;
%         end
%     end
%     xendAPRRtarget=x(end,1)+(y(end,1)-Pnew)/APRR;
%     if y(end,1)<=PtargetEnd
%         yend1=Pnew+dPlower;
%         xend1=xendAPRRtarget+(yend1-Pnew)/APRR;
%         yend2=y(end,1)-dPupper;
%         xend2=xendAPRRtarget-(y(end,1)-yend2)/APRR;
%         if FillTime>=xend2 && FillTime<=xend1
%             PassFillTime='Yes |';
%         else
%             PassFillTime='No  |';
%         end
%     else
%         PassFillTime='No  |';
%     end
% end
% 
% 
% %Tolerance for Leaktime
% tleaknew=0;
% if exist('Ptarget1','var')
%     for i=length(x)-2:-2:2
%         tleak=x(i+1)- x(i);
%         tleaknew=tleaknew+tleak;
%     end
% else
%     for i=length(x)-1:-2:2
%         tleak=(x(i)-(x(i-1)));
%         tleaknew=tleaknew+tleak;
%     end
% end
% 
% LeaktimeTol=30;
% if tleaknew==0
%     PassLeakTime='   NoCheck|';
% else
%     if tleaknew<=LeaktimeTol
%         PassLeakTime='Yes |';
%     else
%         PassLeakTime='No  |';
%     end
% end
% 
% PTolup=PtargetEnd*1.05;
% PToldown=PtargetEnd*0.95;
% 
% if (y(end))>=PToldown && (y(end))<=PTolup
%     PassPtarget='Yes |';
% else
%     PassPtarget='No  |';
% end
% 
% % ind=length(TCU_H2Soc_t0);
% % ind2=round(ind/10);
% % SOC=zeros(ind2,1);
% % for i=ind2:-1:1
% %     SOC(i,1)=TCU_H2Soc_t0(ind-i,1);
% % end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Needs to be updated
% if SOCcalc(end)>=95 && SOCcalc(end)<=100
%     PassSOC='Yes |';
% else
%     PassSOC='No  |';
% end
% SOCtarget='95<x<100';
% 
% 
% %% Show Table in Command Window
% xendAPRRtarget=num2str((xendAPRRtarget-x(1,1))/60);
% FillTime=num2str((FillTime-x(1,1))/60);
% LeaktimeTol=num2str(LeaktimeTol);
% tleaknew=num2str(tleaknew);
% PercentageTol=num2str(PercentageTol);
% Ptarget=num2str(PtargetEnd);
% Pactual=num2str(y(end,1));
% SOC=num2str(SOCcalc(end));
% result=num2str(result);
% 
% fprintf('|                      |Setpoint|      |ActualValue|      |Pass|\n')
% fprintf(['|Filltime+-Tol [Min]   |',xendAPRRtarget,'  |      |',FillTime,'     |      |',PassFillTime,'\n' ] )
% fprintf(['|Leakchecktime [s]     |',LeaktimeTol,'      |      |',tleaknew,'       |      |',PassLeakTime '\n'])
% fprintf(['|Deviations [%%]        |',PercentageTol,'       |      |',result,'         |      ',PassPercentage '\n'])
% fprintf(['|P_Target+-5%% [MPa]    |',Ptarget,' |      |',Pactual,'       |      |',PassPtarget,'\n' ] )
% fprintf(['|SOC [%%]               |',SOCtarget,'|      |',SOC,'         |      |',PassSOC,'\n' ] )
% if exist('resultSFDT','var')
%     fprintf(['|Deviations SFDT [%%]   |',PercentageTol,'       |      |',resultSFDT,'         |      ',PassPercentageSFDT '\n'])
% end
% % if exist('d','var') && exist('xend','var') && xend == x(d+1,1)
% %     fprintf(2, '\n!!! Fill incomplete or wrong User-Input!!!\n')
% % end


end


