%% Combined for all the patient and exercises -------- Need some modifications according to input data

start_path = fullfile('/Users/shehjarsadhu/Desktop/Day2MatlabData');
% Ask user to confirm the folder, or change it.
uiwait(msgbox('Pick a starting folder which has combination of all the files!'));
topLevelFolder = uigetdir(start_path);

filePattern = sprintf('%s/**/*csv', topLevelFolder);
all = dir(filePattern);
j=2;
allFileInfo=all(2);
for i =3:length(all)
    if (all(i).bytes<300) 
        continue;
    else
        allFileInfo(j) = all(i);
        j=j+1;
    end
end

%for normalization
% x = randn(10, 1);
% y = normalize(x, 'range', [0 255]);
%use this code
F=[];
for j=0:6:(length(allFileInfo)-1)
    N = 6;
    Table=[];
    for i = j+1:j+6
        R=[];
        R = readtable(allFileInfo(i).name);
        R=R(:,4:end);
        R=table2array(R);
        %normalization and low pass filtering
        k=size(R,2);
        R=resample(R,1280,size(R,1));
        for z=1:k
            if z == 4 || z ==5 || z == 6 || z == 7 || z == 8 || z == 9
                R(:,z) = normalize(R(:,z), 'range', [0 255]);
            end
            R(:,z) = lowpass(R(:,z),0.4,128); %specify fpass as the centre value
            R(:,z) = (R(:,z)-mean(R(:,z)));  %(R(:,i)-mean(R(:,i)))/std(R(:,i)); - z-score normalization
        end
        %R=array2table(R);
        Table=[Table R];
    end
    Table=array2table(Table);
    %writetable(T,"C:\Users\Apu\Desktop\kaya app\ran\1\combined.csv");
    % Basic combined data file import and normalization and lowpass and then feature extraction file made with features
    % use dlmread(filename,',',1,3)
    % Our sampling frequency=119Hz but take the nearest 2 power which is in this case 128
    data=Table;
    %data = importdata('combined.csv');
    Fs = 128;
    
    ExerciseNames=["Resting Hands On Thigh";"Hands Hold Out";"Finger to Nose";"Finger Tap";"Closed Grip";"Hand Flip"];
    L = 1280; %128*10 because doing each experiment for 10sec 
    T = 1/Fs; %sample period
    t = (0:L-1)*T; %time bins
    n = 2^nextpow2(L); %for zero padding
    f = 0:(Fs/n):(Fs/2-Fs/n); %frequency bins
    dim = 2; %signal along column 

    % Features to be extracted
    % Peak analysis (freq) : perform fft and identify the dominant frequency
    data=table2array(data);
    Y = fft(data',n,dim);
    P2 = abs(Y/L);
    P1 = P2(:,1:n/2+1);
    P1(:,2:end-1) = 2*P1(:,2:end-1);

    for i = 1:size(data,2)
       [pks,locs]=findpeaks(P1(i,:));
       [M,I]=max(pks);
       DomFreqMag(i,1) = P1(i,locs(I));
       DomFreq(i,1) = f(locs(I));
       Label{i,1} = 'Exercises';
    end
%     z=1;
%     figure
%     for i=1:60
%        subplot(6,10,i)
%        plot(f,P1(i,1:n/2))
%        xlim([0,10])
%        title(['ExNo.',num2str(z+floor((i-1)/10))])
%     end
%     sgtitle('Test Person No.%d ~FFT plots for all (Thumb,Ring,Middle,Index,Accx,Accy,Accz,Gyrx,Gyry,Gyrz) in specific order',ceil((j+1)/6))
    % spectral measurements : signal power, bandwidth, mean frequency, median frequency
    for i = 1:size(data,2)
       PWR(i,1) = bandpower(data(:,i),Fs,[0 10]);
       OBW(i,1) = obw(data(:,i),Fs);
       MeanFreq(i,1) = meanfreq(data(:,i),Fs);
       MedFreq(i,1) = medfreq(data(:,i),Fs);
       MaxValue(i,1) = max(data(:,i));
       MinValue(i,1) = min(data(:,i));
       Var(i,1) = var(data(:,i));
    end

    % Peak analysis (time) : max peak, min peak, mean peak, variance in peaks
    for i = 1:size(data,2)
       [pks,locs]=findpeaks(data(i,:));
       [M,I]=max(pks);
       MaxPeak(i,1) = data(i,locs(I));
       [M,I]=min(pks);
       MinPeak(i,1) = data(i,locs(I));
       MeanPeak(i,1) = mean(pks);
       VarPeak(i,1) = var(pks);
    end

    %table for all the features variance peak.
    for i=1:9:size(data,2)
        VarPThumb(floor(i/9)+1,1)=VarPeak(i,1);
        VarPIndex(floor(i/9)+1,1)=VarPeak(i+1,1);
        %VarPMiddle(floor(i/10)+1,1)=VarPeak(i+2,1);
        VarPRing(floor(i/9)+1,1)=VarPeak(i+2,1);
        VarPAccx(floor(i/9)+1,1)=VarPeak(i+3,1);
        VarPAccy(floor(i/9)+1,1)=VarPeak(i+4,1);
        VarPAccz(floor(i/9)+1,1)=VarPeak(i+5,1);
        VarPGyrx(floor(i/9)+1,1)=VarPeak(i+6,1);
        VarPGyry(floor(i/9)+1,1)=VarPeak(i+7,1);
        VarPGyrz(floor(i/9)+1,1)=VarPeak(i+8,1);
    end
    
    TVarP=table(VarPThumb, VarPIndex, VarPRing, VarPAccx, VarPAccy, VarPAccz, VarPGyrx, VarPGyry, VarPGyrz);
    
    for i=1:9:size(data,2)
        MeanPThumb(floor(i/9)+1,1)=MeanPeak(i,1);
        MeanPIndex(floor(i/9)+1,1)=MeanPeak(i+1,1);
        %MeanPMiddle(floor(i/10)+1,1)=MeanPeak(i+2,1);
        MeanPRing(floor(i/9)+1,1)=MeanPeak(i+2,1);
        MeanPAccx(floor(i/9)+1,1)=MeanPeak(i+3,1);
        MeanPAccy(floor(i/9)+1,1)=MeanPeak(i+4,1);
        MeanPAccz(floor(i/9)+1,1)=MeanPeak(i+5,1);
        MeanPGyrx(floor(i/9)+1,1)=MeanPeak(i+6,1);
        MeanPGyry(floor(i/9)+1,1)=MeanPeak(i+7,1);
        MeanPGyrz(floor(i/9)+1,1)=MeanPeak(i+8,1);
    end
    TMeanP=table(MeanPThumb, MeanPIndex, MeanPRing, MeanPAccx, MeanPAccy, MeanPAccz, MeanPGyrx, MeanPGyry, MeanPGyrz);
    for i=1:9:size(data,2)
        MinPThumb(floor(i/9)+1,1)=MinPeak(i,1);
        MinPIndex(floor(i/9)+1,1)=MinPeak(i+1,1);
        %MinPMiddle(floor(i/9)+1,1)=MinPeak(i+2,1);
        MinPRing(floor(i/9)+1,1)=MinPeak(i+2,1);
        MinPAccx(floor(i/9)+1,1)=MinPeak(i+3,1);
        MinPAccy(floor(i/9)+1,1)=MinPeak(i+4,1);
        MinPAccz(floor(i/9)+1,1)=MinPeak(i+5,1);
        MinPGyrx(floor(i/9)+1,1)=MinPeak(i+6,1);
        MinPGyry(floor(i/9)+1,1)=MinPeak(i+7,1);
        MinPGyrz(floor(i/9)+1,1)=MinPeak(i+8,1);
    end
    TMinP=table(MinPThumb, MinPIndex, MinPRing, MinPAccx, MinPAccy, MinPAccz, MinPGyrx, MinPGyry, MinPGyrz);
    for i=1:9:size(data,2)
        MaxPThumb(floor(i/9)+1,1)=MaxPeak(i,1);
        MaxPIndex(floor(i/9)+1,1)=MaxPeak(i+1,1);
        %MaxPMiddle(floor(i/9)+1,1)=MaxPeak(i+2,1);
        MaxPRing(floor(i/9)+1,1)=MaxPeak(i+2,1);
        MaxPAccx(floor(i/9)+1,1)=MaxPeak(i+3,1);
        MaxPAccy(floor(i/9)+1,1)=MaxPeak(i+4,1);
        MaxPAccz(floor(i/9)+1,1)=MaxPeak(i+5,1);
        MaxPGyrx(floor(i/9)+1,1)=MaxPeak(i+6,1);
        MaxPGyry(floor(i/9)+1,1)=MaxPeak(i+7,1);
        MaxPGyrz(floor(i/9)+1,1)=MaxPeak(i+8,1);
    end
    TMaxP=table(MaxPThumb, MaxPIndex, MaxPRing, MaxPAccx, MaxPAccy, MaxPAccz, MaxPGyrx, MaxPGyry, MaxPGyrz);
    for i=1:9:size(data,2)
        DFMThumb(floor(i/9)+1,1)=DomFreqMag(i,1);
        DFMIndex(floor(i/9)+1,1)=DomFreqMag(i+1,1);
        %DFMMiddle(floor(i/9)+1,1)=DomFreqMag(i+2,1);
        DFMRing(floor(i/9)+1,1)=DomFreqMag(i+2,1);
        DFMAccx(floor(i/9)+1,1)=DomFreqMag(i+3,1);
        DFMAccy(floor(i/9)+1,1)=DomFreqMag(i+4,1);
        DFMAccz(floor(i/9)+1,1)=DomFreqMag(i+5,1);
        DFMGyrx(floor(i/9)+1,1)=DomFreqMag(i+6,1);
        DFMGyry(floor(i/9)+1,1)=DomFreqMag(i+7,1);
        DFMGyrz(floor(i/9)+1,1)=DomFreqMag(i+8,1);
    end
    TDFM=table(ExerciseNames, DFMThumb, DFMIndex, DFMRing, DFMAccx, DFMAccy, DFMAccz, DFMGyrx, DFMGyry, DFMGyrz);
    for i=1:9:size(data,2)
        DFThumb(floor(i/9)+1,1)=DomFreq(i,1);
        DFIndex(floor(i/9)+1,1)=DomFreq(i+1,1);
        %DFMiddle(floor(i/9)+1,1)=DomFreq(i+2,1);
        DFRing(floor(i/9)+1,1)=DomFreq(i+2,1);
        DFAccx(floor(i/9)+1,1)=DomFreq(i+3,1);
        DFAccy(floor(i/9)+1,1)=DomFreq(i+4,1);
        DFAccz(floor(i/9)+1,1)=DomFreq(i+5,1);
        DFGyrx(floor(i/9)+1,1)=DomFreq(i+6,1);
        DFGyry(floor(i/9)+1,1)=DomFreq(i+7,1);
        DFGyrz(floor(i/9)+1,1)=DomFreq(i+8,1);
    end
    TDF=table(DFThumb, DFIndex, DFRing, DFAccx, DFAccy, DFAccz, DFGyrx, DFGyry, DFGyrz);
    for i=1:9:size(data,2)
        PThumb(floor(i/9)+1,1)=PWR(i,1);
        PIndex(floor(i/9)+1,1)=PWR(i+1,1);
        %PMiddle(floor(i/9)+1,1)=PWR(i+2,1);
        PRing(floor(i/9)+1,1)=PWR(i+2,1);
        PAccx(floor(i/9)+1,1)=PWR(i+3,1);
        PAccy(floor(i/9)+1,1)=PWR(i+4,1);
        PAccz(floor(i/9)+1,1)=PWR(i+5,1);
        PGyrx(floor(i/9)+1,1)=PWR(i+6,1);
        PGyry(floor(i/9)+1,1)=PWR(i+7,1);
        PGyrz(floor(i/9)+1,1)=PWR(i+8,1);
    end
    TP=table(PThumb, PIndex, PRing, PAccx, PAccy, PAccz, PGyrx, PGyry, PGyrz);
    for i=1:9:size(data,2)
        MeanFThumb(floor(i/9)+1,1)=MeanFreq(i,1);
        MeanFIndex(floor(i/9)+1,1)=MeanFreq(i+1,1);
        %MeanFMiddle(floor(i/9)+1,1)=MeanFreq(i+2,1);
        MeanFRing(floor(i/9)+1,1)=MeanFreq(i+2,1);
        MeanFAccx(floor(i/9)+1,1)=MeanFreq(i+3,1);
        MeanFAccy(floor(i/9)+1,1)=MeanFreq(i+4,1);
        MeanFAccz(floor(i/9)+1,1)=MeanFreq(i+5,1);
        MeanFGyrx(floor(i/9)+1,1)=MeanFreq(i+6,1);
        MeanFGyry(floor(i/9)+1,1)=MeanFreq(i+7,1);
        MeanFGyrz(floor(i/9)+1,1)=MeanFreq(i+8,1);
    end
    TMeanF=table(MeanFThumb, MeanFIndex, MeanFRing, MeanFAccx, MeanFAccy, MeanFAccz, MeanFGyrx, MeanFGyry, MeanFGyrz);
    for i=1:9:size(data,2)
        MedianFThumb(floor(i/9)+1,1)=MedFreq(i,1);
        MedianFIndex(floor(i/9)+1,1)=MedFreq(i+1,1);
        %MedianFMiddle(floor(i/9)+1,1)=MedFreq(i+2,1);
        MedianFRing(floor(i/9)+1,1)=MedFreq(i+2,1);
        MedianFAccx(floor(i/9)+1,1)=MedFreq(i+3,1);
        MedianFAccy(floor(i/9)+1,1)=MedFreq(i+4,1);
        MedianFAccz(floor(i/9)+1,1)=MedFreq(i+5,1);
        MedianFGyrx(floor(i/9)+1,1)=MedFreq(i+6,1);
        MedianFGyry(floor(i/9)+1,1)=MedFreq(i+7,1);
        MedianFGyrz(floor(i/9)+1,1)=MedFreq(i+8,1);
    end
    TMedianF=table(MedianFThumb, MedianFIndex, MedianFRing, MedianFAccx, MedianFAccy, MedianFAccz, MedianFGyrx, MedianFGyry, MedianFGyrz);
    for i=1:9:size(data,2)
        BWThumb(floor(i/9)+1,1)=OBW(i,1);
        BWIndex(floor(i/9)+1,1)=OBW(i+1,1);
        %BWMiddle(floor(i/9)+1,1)=OBW(i+2,1);
        BWRing(floor(i/9)+1,1)=OBW(i+2,1);
        BWAccx(floor(i/9)+1,1)=OBW(i+3,1);
        BWAccy(floor(i/9)+1,1)=OBW(i+4,1);
        BWAccz(floor(i/9)+1,1)=OBW(i+5,1);
        BWGyrx(floor(i/9)+1,1)=OBW(i+6,1);
        BWGyry(floor(i/9)+1,1)=OBW(i+7,1);
        BWGyrz(floor(i/9)+1,1)=OBW(i+8,1);
    end
    TBW=table(BWThumb, BWIndex, BWRing, BWAccx, BWAccy, BWAccz, BWGyrx, BWGyry, BWGyrz);
    for i=1:9:size(data,2)
        MaxValueThumb(floor(i/9)+1,1)=MaxValue(i,1);
        MaxValueIndex(floor(i/9)+1,1)=MaxValue(i+1,1);
        %MaxValueMiddle(floor(i/9)+1,1)=MaxValue(i+2,1);
        MaxValueRing(floor(i/9)+1,1)=MaxValue(i+2,1);
        MaxValueAccx(floor(i/9)+1,1)=MaxValue(i+3,1);
        MaxValueAccy(floor(i/9)+1,1)=MaxValue(i+4,1);
        MaxValueAccz(floor(i/9)+1,1)=MaxValue(i+5,1);
        MaxValueGyrx(floor(i/9)+1,1)=MaxValue(i+6,1);
        MaxValueGyry(floor(i/9)+1,1)=MaxValue(i+7,1);
        MaxValueGyrz(floor(i/9)+1,1)=MaxValue(i+8,1);
    end
    TMaxValue=table(MaxValueThumb, MaxValueIndex, MaxValueRing, MaxValueAccx, MaxValueAccy, MaxValueAccz, MaxValueGyrx, MaxValueGyry, MaxValueGyrz);
    for i=1:9:size(data,2)
        MinValueThumb(floor(i/9)+1,1)=MinValue(i,1);
        MinValueIndex(floor(i/9)+1,1)=MinValue(i+1,1);
        %MinValueMiddle(floor(i/9)+1,1)=MinValue(i+2,1);
        MinValueRing(floor(i/9)+1,1)=MinValue(i+2,1);
        MinValueAccx(floor(i/9)+1,1)=MinValue(i+3,1);
        MinValueAccy(floor(i/9)+1,1)=MinValue(i+4,1);
        MinValueAccz(floor(i/9)+1,1)=MinValue(i+5,1);
        MinValueGyrx(floor(i/9)+1,1)=MinValue(i+6,1);
        MinValueGyry(floor(i/9)+1,1)=MinValue(i+7,1);
        MinValueGyrz(floor(i/9)+1,1)=MinValue(i+8,1);
    end
    TMinValue=table(MinValueThumb, MinValueIndex, MinValueRing, MinValueAccx, MinValueAccy, MinValueAccz, MinValueGyrx, MinValueGyry, MinValueGyrz);
    for i=1:9:size(data,2)
        VarThumb(floor(i/9)+1,1)=Var(i,1);
        VarIndex(floor(i/9)+1,1)=Var(i+1,1);
        %VarMiddle(floor(i/9)+1,1)=Var(i+2,1);
        VarRing(floor(i/9)+1,1)=Var(i+2,1);
        VarAccx(floor(i/9)+1,1)=Var(i+3,1);
        VarAccy(floor(i/9)+1,1)=Var(i+4,1);
        VarAccz(floor(i/9)+1,1)=Var(i+5,1);
        VarGyrx(floor(i/9)+1,1)=Var(i+6,1);
        VarGyry(floor(i/9)+1,1)=Var(i+7,1);
        VarGyrz(floor(i/9)+1,1)=Var(i+8,1);
    end
    TVar=table(VarThumb, VarIndex, VarRing, VarAccx, VarAccy, VarAccz, VarGyrx, VarGyry, VarGyrz);
    
    %appending all the tables together
    finaltable = [TDFM TDF TP TMeanF TMedianF TBW TMaxP TMinP TMeanP TVarP TMaxValue TMinValue TVar];
    %writetable(finaltable,"AllFeatureTable.csv");
    %finaltable = table2array(finaltable);
    F=[F ; finaltable];
end
%F=array2table(F);
writetable(F,"/Users/shehjarsadhu/Desktop/AllfeaturesDay2_onlyright.csv");