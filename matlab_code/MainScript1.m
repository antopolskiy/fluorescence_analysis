# 1. Clustering for the heatmap order of cells
# 2. Clip at zero when calculating the AUC
# 3. Maximum amplitude of of the response to each stimulus
# 4. Table of responses to the stimuli


%% Load xlsx file into the matlab variables as a numerical matrix
##dF = readtable('S18.xlsx'); manually open the xlsx file of interest and
% sellect the range excluding frame number. assign it to the variable dF
%% Define path to the custom functions
##path(path,genpath('/Users/madinoid/Desktop/MadinaMatLab/MatLab functions'));

%df = load('dF_Variable.mat');

%% Convert Frames into time(min), define variables frameNum,cellNum,frames,time
frameNum = size(dF,1);
cellNum = size(dF,2);
frames = (1:frameNum)';
##prompt = 'input the frame duration t in sec';
##frameDur = input(prompt);
frameDur = 361
time = frameDur*frames/60';
Cell=zeros(cellNum,1);
for i=1:cellNum
    Cell(i,1)=i;
end
%% Convert row traces into dF/F
##prompt = 'how many frames is your baseline?';
##F0N = input(prompt);
F0N = 161
F0 = repmat(mean(dF(1:F0N,:),1),frameNum,1);

%horizontal vector where each cell is the average of the first F0N rows of the corresponding dF matrix columns
dF_F = (dF-F0)./F0*100;
% dF_F_raw = dF_F;
%% draw vertical line indicating where the stimulus was applied, 
% prompt user's input
##prompt = 'type how many stimuli were tested in this trial?...';
##stimNum = input(prompt);
stimNum = 3
s_Array = zeros(stimNum,3);

s_on = [162, 307, 401]
for i = 1 : stimNum
##query=['what frame is onset of stim number ',num2str(i)];
##disp(query);
##prompt = '? ';
##s_on = input(prompt);
% query=['what frame is endpoint for stim number (input 0 if none) ',num2str(i)];
% disp(query);
% prompt = '? ';
% s_off = input(prompt);
s_Array(i,1)=i;
s_Array(i,2)=s_on(i);
% s_Array(i,3)=s_off;
end

%% Visualize 1 trace at a time for the QC
% This plots verticle lines for every x value defined by the user. for each stimulus
f=figure(1);
h=plot(time,dF_F);
s_Array_min = s_Array(1:end,2:3)*frameDur/60;
s_onLine=vline(s_Array_min(:,1));
% s_offLine=vline(s_Array_min(:,2));

% This plots one trace per figure for QC 
handles=struct('f',f,'h',h);
handles.visible_index=1;%keep track of which plot is selected
set(h,'Visible','off');
set(h(handles.visible_index),'Visible','on');
guidata(f,handles)
set(f,'KeyPressFcn',@SetVisibility)
%% Plot heatmap of dF/F traces over time
f2 = figure;
fh = heatmap((dF_F)','XLabel','Time(min)','YLabel','Cell #','Colormap',jet);

%% the following calculates the mean and stand dev of the baseline and generates 2 arrays,
% fbase ("baseline fluorescence") and fstdev (stand deviation of baseline fluorescence)
fbase = mean(dF_F(1:F0N,:));
fstdev = std(dF_F(1:F0N,:));

%% the following section looks for any points in the baseline (1st F0N points) that are > (mean
%baseline + 2x stand dev) and replaces that value with the mean baseline fluorescence. 
% That is, this removes any large signals in the baseline (1st F0N points)
for i=1:cellNum
  for row=1:F0N
    if(dF_F(row,i)>(fbase(1,i))+2*fstdev(1,i))
      dF_F(row,i)=fbase(1,i);
    end
  end
end
f3 = figure;
fh = heatmap((dF_F)','XLabel','Time(min)','YLabel','Cell #','Colormap',jet);

%% the following recalculates a new mean and new stand deviation for the
%baseline (q points), having removed any large spontaneous signals from
%the baseline. Then this calculates a "criterion" for each cell that is
%equal to the baseline of that cell + m standard deviations. The user has
%input m at the beginning
newfbase = mean(dF_F(1:frameNum,:));
newfstdev = std(dF_F(1:frameNum,:));
criterion=newfbase+2*newfstdev;
%% NEW NEEDS TO BE CHECKED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

data = transpose(dF_F);
multiplier1 = 1.25;
multiplier2 = 3;

NoDriftdata=zeros(cellNum,frameNum); %create a matrix of zeros with the number of cells = #rows, and number of columns = #scan
for x=1:cellNum
    for c=1:20
        NoDriftdata(x,c)=data(x,c);
    end
    factor=multiplier1*std(data(x,1:30)); %this resets "factor" to be a multiplicand of the st dev of baseline
    for c=21:frameNum-1
        currentpoint=data(x,c);
        previouspoint=data(x,c-1);
        previousnodrift=NoDriftdata(x,c-1);
        abs1=abs(currentpoint-previouspoint);
        abs2=abs(currentpoint-previousnodrift);
        if abs1 > factor || abs2 > factor  %this says if  abs1 > factor OR abs2 > factor then execute the next line
            NoDriftdata(x,c)=mean(NoDriftdata(x,(c-4):(c-1))); %creates a matrix "NoDrift" that is dF_F w/o peaks (i.e. has eliminated stim-evoked responses)
        else  NoDriftdata(x,c)=data(x,c);
        end
    end
    movingav(x,:)=movmean(NoDriftdata(x,:),10);
    result(x,:)=data(x,:)-movingav(x,:); %creates a matrix, "result", that has all the drift-corrected data.  Result is the equivalent of dF_F (raw data) or "data" (also, raw data)
    transresult1=transpose(result);

    factor=multiplier2*std(data(x,1:30)); %this resets "factor" to be a multiplicand of the st dev of baseline
    for c=21:frameNum-1
        currentpoint=data(x,c);
        previouspoint=data(x,c-1);
        previousnodrift=NoDriftdata(x,c-1);
        abs1=abs(currentpoint-previouspoint);
        abs2=abs(currentpoint-previousnodrift);
        if abs1 > factor || abs2 > factor  %this says if  abs1 > factor OR abs2 > factor then execute the next line
            NoDriftdata(x,c)=mean(NoDriftdata(x,(c-4):(c-1))); %creates a matrix "NoDrift" that is dF_F w/o peaks (i.e. has eliminated stim-evoked responses)
        else  NoDriftdata(x,c)=data(x,c);
        end
    end
    movingav(x,:)=movmean(NoDriftdata(x,:),10);
    result(x,:)=data(x,:)-movingav(x,:); %creates a matrix, "result", that has all the drift-corrected data.  Result is the equivalent of dF_F (raw data) or "data" (also, raw data)

transresult2=transpose(result);

end
transresult3=transpose(result);

%This section plots each cell 
f4=figure(4);
y2=transresult1+ 300;
y3=transresult2+ 600;
y4=transresult3+ 900;
h4=plot(time,dF_F,'b',time,y2,'k',time,y3,'r');
h4=reshape(h4,cellNum,3);
s_Array_min = s_Array(1:end,2:3)*frameDur/60;
s_onLine=vline(s_Array_min(:,1));
handles=struct('f',f4,'h',h4);
handles.visible_index=1;
set(h4,'Visible','off');
set(h4(handles.visible_index,:),'Visible','on');
guidata(f4,handles)
set(f4,'KeyPressFcn',@SetVisibility)

%%
%This section plots each cell with blue=original data,black=drift,  red=drift eliminated
f4=figure(4);
y2=movingav'+ 300;
y3=result'+ 600;
h4=plot(time,dF_F,'b',time,y2,'k',time,y3,'r');
h4=reshape(h4,cellNum,3);
s_Array_min = s_Array(1:end,2:3)*frameDur/60;
s_onLine=vline(s_Array_min(:,1));
handles=struct('f',f4,'h',h4);
handles.visible_index=1;
set(h4,'Visible','off');
set(h4(handles.visible_index,:),'Visible','on');
guidata(f4,handles)
set(f4,'KeyPressFcn',@SetVisibility)

%% SELECT FACTOR TO REDUCE DRIFT
% restart from HERE if you need to abort the plotting (hit CNTRL/C, or close fig window)
% the following section attempts to eliminate drift or "wobble" in the baseline
% this creates a matrix, "result" that has all the drift-corrected values
data = transpose(dF_F);
query='what is your cutoff factor for reducing drift and wobble (i.e. factor x st dev of baseline),e.g. try 1.25 to 3';
disp(query);
prompt = '? ';
multiplier = input(prompt);
NoDriftdata=zeros(cellNum,frameNum); %create a matrix of zeros with the number of cells = #rows, and number of columns = #scan
for x=1:cellNum
    for c=1:20
        NoDriftdata(x,c)=data(x,c);
    end
    factor=multiplier*std(data(x,1:30)); %this resets "factor" to be a multiplicand of the st dev of baseline
    for c=21:frameNum-1
        currentpoint=data(x,c);
        previouspoint=data(x,c-1);
        previousnodrift=NoDriftdata(x,c-1);
        abs1=abs(currentpoint-previouspoint);
        abs2=abs(currentpoint-previousnodrift);
        if abs1 > factor || abs2 > factor  %this says if  abs1 > factor OR abs2 > factor then execute the next line
            NoDriftdata(x,c)=mean(NoDriftdata(x,(c-4):(c-1))); %creates a matrix "NoDrift" that is dF_F w/o peaks (i.e. has eliminated stim-evoked responses)
        else  NoDriftdata(x,c)=data(x,c);
        end
    end
    movingav(x,:)=movmean(NoDriftdata(x,:),10);
    result(x,:)=data(x,:)-movingav(x,:); %creates a matrix, "result", that has all the drift-corrected data.  Result is the equivalent of dF_F (raw data) or "data" (also, raw data)
end
transresult=transpose(result);
%% the following calculates the area under the curve (AUC) for uncorrected 
% traces, between each stimulus onset to the point defined by the user.
% column1 = cell number
% column2 = baseline AUC
% columns2:stimNum = stimulus-specific AUC 
AUC=zeros([cellNum stimNum+2]);
prompt = 'How many frames after stimulus onset do you want to consider for AUC calculation?...';
framesAUC = input(prompt);
for i=1:cellNum
    for j=1:stimNum
        begin=s_Array(j,2);
        final=begin + framesAUC;
        AUC(i,1)=i;            %this numbers all the cells
        AUC(i,2)=trapz(dF_F(1:framesAUC,i)); %this calcultaes baseline AUC (frame1 - #frames defined by the user)
        AUC(i,j+2)=trapz(dF_F(begin:final,i));%this calculates AUC for each cell following each stimulus
        
    end
end


%% the following does the same as above, but for the traces that have been corrected for drift
AUCnodrift=zeros([cellNum stimNum+2]);
for i=1:cellNum
    for j=1:stimNum
        begin=s_Array(j,2);
        final=begin + framesAUC;
        AUCnodrift(i,1)=i;
        AUC(i,2) = trapz(transresult(1:framesAUC,i));
        AUCnodrift(i,j+1)=trapz(transresult(begin:final,i));
    end
end

% %% the following finds the peak response after each stimulus for the uncorrected traces
% MAX=zeros(cellNum, stimNum+1);
% for i=1:cellNum
%         MAX(i,1)=i;
%         for j=1:stimNum
%         begin=s_Array(j,2);
%         final=begin + framesAUC;
%         MAX(i,j+1)=max(dF_F(begin:final,i));
%         end
% end
% 
% %% the following finds the peak response after each stimulus for the drift-corrected traces
% for i=1:cellNum
%     for j=1:stimNum
%         begin=s_Array(j,2);
%         final=begin + framesAUC;
%         MAXnodrift(i,1)=i;
%         MAXnodrift(i,j+1)=max(transresult(begin:final,i));
%     end
% end
    
%% the following creates a new matrix (for uncorrected traces): (AUC)trace > (AUC)baseline
% 1 - indicates that the cell responded to the stimulus
% 0 - indicates no responce
% treshold(t) is defined by the user as the number of SD above the baseline

prompt = 'How many SD above the baseline is the responce treshold?...';
t_Number = input(prompt);
TH = abs(t_Number * newfstdev);
responders = (AUC(:,2)+TH') < AUC(:,setdiff(1:size(AUC,2),[1 2]));
[respNum, stim] = find(responders);%% this extracts indeces of the responding cells (2 vectors: rows and columns) 
respNum_unique = unique(respNum);
respAUC = AUC(respNum_unique,2:end);
figure
bp = boxplot(respAUC);
%% TO DO: this extracts an AUC for the traces of the correspondinc indeces
% NEED TO GENERATE A 3D ARRAY EXTRACTING TRACES OF THE RESPONDING CELLS FOR
% ALL STIMULI. 1D = frmaeNum, 2D = cellNum, 3D = stimNum


%% the following creates a new matrix (for corrected traces): (AUC)trace > (AUC)baseline
% 1 - indicates that the cell responded to the stimulus
% 0 - indicates no responce
% treshold(t) is defined by the user as the number of SD above the baseline

responders_nodrift = (AUCnodrift(:,2)+TH') < AUCnodrift(:,setdiff(1:size(AUCnodrift,2),[1 2]));
[respNum1, stim1] = find(responders_nodrift);
respNum1_unique = unique(respNum1);
respAUC1 = AUCnodrift(respNum1_unique,2:end);
figure
bp = boxplot(respAUC1);
%% the following generates an output cell array containing all the measured parameters 
for j = 1:stimNum
    query=['what is the name of stim number followed by the date of the experiment ',num2str(j)];
    disp(query);
    prompt = '? ';
    val1 = input(prompt,'s'); %stimulus name
    val2 = sum(responders(:,j)); % number of cells responding to the stimulus
    val3 = sum(responders_nodrift(:,j)); % number of cells responding to the stimulus (drift corrected)
    val4 = cellNum; %total number of cells analyzed in this experiment
    c_Array(j,:) = {val1, val2, val3, val4};
end
%% Saving variables
prompt = 'Enter the date of the experiment or how you want to name this file...';
eName = input(prompt,'s');
save(eName)
%
% %% This script is written for manually sellected data, it uses variables 
% % generated by the MainScript. Only run it for the same series that were
% % analyzed by the MainScript.
% % 
% % This code will:
% %    1) Plot staked traces: figure showing traces of all the cells
% %    responding to this stimulus;
% %    3) Count the number of cells responding to this stimulus, % of reponding cells;
% %    2) Calculate AUC for the stimulus vs baseline for the number of frames
% %    sellected by the user.
% for k = 1:stimNum    
% %% Define path to the master table
% load('master');
% %% Import xlsx file with the dF_F values for the cells of interest for each
% % stimulus. Name the file SdF_F.
% %% This plots the data along with v-lines when the stimulus was applied
% prompt = 'Enter the date of the experiment followed the stimulus  that is being analyzed?...';
% sName = input(prompt,'s');
% SF = figure('NumberTitle', 'off', 'Name', sName);
% SFh = stackedplot(time,SdF_F);
% s_onLine=vline(s_Array_min(:,1));
% s_offLine=vline(s_Array_min(:,2));
% 
% %% This generates a horizonatl array with stimulus name, number of responded 
% % cells, total number of cells, % responded cells.
% val1 = sName;
% val2 = size(SdF_F,2);
% val3 = cellNum;
% val4 = size(SdF_F,2)/cellNum*100;
% S_countArray = {val1, val2, val3, val4};
% %% this will append S_countArray into a master array table that has the number
% % of responding cells from multiple experiments. Saves all the variables
% master = [master; S_countArray];
% save(sName);
% save(sName,'master',version);           
% end
%% Split traces into categories based on AUC for each stimulus vs baseline AUC
%f3 = figure(3);
%sp = stackedplot(time,dF_F); 
%onLine=vline(s_Array_min(:,1));
%for i = 1 : stimNum
%prompt = ['Enter cells that respond to stimulus ',num2str(i),'[n1 n2 ...ni] ADD 0s TO MAKE IMPUTS EVEN IN LENGTH...'];
%sv{i} = input(prompt); % string of inputs for each stimulus
%end
%sm = cell2mat(sv');% converts each variable of sv_dialog into a separate vector
%% plot traces for each stimulus based on the visualized QC 
%figure (4);
%for i = 1:stimNum;
   % subplot (2,2,(2/i));
   % stackedplot(time,dF_F(:,sm(i,:)));
%end