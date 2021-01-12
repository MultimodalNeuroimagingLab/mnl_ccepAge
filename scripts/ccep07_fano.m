%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end

%% connections from one region to another

region_start = input('Choose roi where connections start [temporal, frontal, parietal, occipital]: ','s');
region_end = input('Choose roi where connections end [temporal, frontal, parietal, occipital]: ','s');

% region_start = 'occipital';
% region_end = 'occipital';

out = ccep_connectRegions(n1Latencies,region_start,region_end);

% initialize output: age, mean and variance in latency per subject
my_output = NaN(length(out.sub),3);

% get variable per subject
for kk = 1:length(out.sub)
    my_output(kk,1) = out.sub(kk).age;
    my_output(kk,2) = mean(out.sub(kk).latencies,'omitnan');
    my_output(kk,3) = var(out.sub(kk).latencies,'omitnan');
end

figure
plot(my_output(:,1),1000*my_output(:,2),'k.','MarkerSize',10)
xlabel('age (years)'),ylabel('mean dT (ms)')
[r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
title(['r=' num2str(r,3) ' p=' num2str(p,3)])

%% figure latency connections from one region to another
clear out
regions_connect = {'temporal','central','parietal','frontal'};
conn_matrix = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];

for outInd = 1:size(conn_matrix,1)
    region_start = regions_connect{conn_matrix(outInd,1)};
    region_end = regions_connect{conn_matrix(outInd,2)};
    temp = ccep_connectRegions(n1Latencies,region_start,region_end);
    out{outInd} = temp;
    out{outInd}.name = {region_start,region_end};
end

%% mean vs var


figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),3);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = mean(1000*out{outInd}.sub(kk).latencies,'omitnan');
        my_output(kk,3) = var(1000*out{outInd}.sub(kk).latencies,'omitnan');
    end

    % age vs mean CCEP
    subplot(4,4,outInd),hold on
    plot([1:60],[1:60],'b')
    plot(my_output(:,2),my_output(:,3),'k.','MarkerSize',10)
    xlabel('mean dT (ms)'),ylabel('var dT')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),2),my_output(~isnan(my_output(:,2)),3),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
%     xlim([0 60])%, ylim([0 100])
    
    % Yeatman et al., fit a second order polynomial:
    % y  = w1* age^2 * w2*age + w3
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),2),my_output(~isnan(my_output(:,2)),3),1);
    x_mean = [1:1:60];
    y_fit = P(1)*x_mean + P(2);
    plot(x_mean,y_fit,'r')

end

%%
%%
%% age vs var

figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),3);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = mean(out{outInd}.sub(kk).latencies,'omitnan');
        my_output(kk,3) = var(out{outInd}.sub(kk).latencies,'omitnan');
    end

    % age vs mean CCEP
    subplot(4,4,outInd),hold on
    plot(my_output(:,1),1000*my_output(:,3),'k.','MarkerSize',10)
    xlabel('age (years)'),ylabel('std dT')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),3),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
    xlim([0 60])%, ylim([0 100])
    
    % Let's fit a first order polynomial:
    % y  =  w1*age + w2
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),1),1000*my_output(~isnan(my_output(:,2)),3),1);
    x_age = [1:1:50];
    y_fit = P(1)*x_age + P(2);
    plot(x_age,y_fit,'r')
end

%% age vs fano

figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),4);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = mean(out{outInd}.sub(kk).latencies,'omitnan');
        my_output(kk,3) = var(out{outInd}.sub(kk).latencies,'omitnan');
        my_output(kk,4) = my_output(kk,3)/my_output(kk,2); % fano
    end

    % age vs mean CCEP
    subplot(4,4,outInd),hold on
    plot(my_output(:,1),my_output(:,4),'k.','MarkerSize',10)
    xlabel('age'),ylabel('fano')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),4),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
%     xlim([0 60])%, ylim([0 100])
    
    % Yeatman et al., fit a second order polynomial:
    % y  = w1* age^2 * w2*age + w3
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),4),1);
    x_mean = [1:1:50];
    y_fit = P(1)*x_mean + P(2);
    plot(x_mean,y_fit,'r')

end
