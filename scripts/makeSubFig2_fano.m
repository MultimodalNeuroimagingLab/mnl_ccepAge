%% makeSubFig2_fano
% this code looks into:
% 1. the correlation between mean and variance in latency
% --> a higher mean was correlated to a higher variance in 14 connections
% (p<0.05, FDR corrected)
% 2. the correlation between age and variance in latency
% --> this does not show any significant results 
% 3. the correlation between age and fano factor: variance/mean.
% --> this does not show any significant results



%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end

%% connections from one region to another

region_start = input('Choose roi where connections start [temporal, frontal, parietal, central]: ','s');
region_end = input('Choose roi where connections end [temporal, frontal, parietal, central]: ','s');

% categorize anatomical regions
ccep_categorizeAnatomicalRegions % --> gives roi_name with order of regions

roi_start = find(strcmp(roi_name,region_start));
roi_end = find(strcmp(roi_name,region_end));

out = ccep_connectRegions(n1Latencies,roi{roi_start},roi{roi_end});

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

% categorize anatomical regions
ccep_categorizeAnatomicalRegions % --> gives roi_name with order of regions

regions_connect = {'temporal','central','parietal','frontal'};
conn_matrix = [1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4];

out = cell(1,size(conn_matrix,1));
for outInd = 1:size(conn_matrix,1)
    
    % print what is loaded
    region_start = roi_name{conn_matrix(outInd,1)};
    region_end = roi_name{conn_matrix(outInd,2)};
    string_start = [repmat('%d, ',1,size(roi{conn_matrix(outInd,1)},2)-1), '%d'];
    string_end = [repmat('%d, ',1,size(roi{conn_matrix(outInd,2)},2)-1), '%d'];
    fprintf(['Run %s (' string_start ') - %s (' string_end ') \n'],...
        region_start,str2double(roi{conn_matrix(outInd,1)}),...
        region_end,str2double(roi{conn_matrix(outInd,2)}))
    
    temp = ccep_connectRegions(n1Latencies,roi{conn_matrix(outInd,1)},roi{conn_matrix(outInd,2)});
    out{outInd} = temp;
    out{outInd}.name = {region_start,region_end};

end

%% mean vs var
% calculate correlation and p
p_all = zeros(size(conn_matrix,1),1);
r_all = zeros(size(conn_matrix,1),1);

figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),3);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        
        if length(out{outInd}.sub(kk).latencies) > 1
            my_output(kk,2) = mean(1000*out{outInd}.sub(kk).latencies,'omitnan');
            my_output(kk,3) = var(1000*out{outInd}.sub(kk).latencies,'omitnan');
        end
    end
    
    % age vs mean CCEP
    subplot(4,4,outInd),hold on
%     plot(1:100,1:100,'b')
    plot(my_output(:,2),my_output(:,3),'k.','MarkerSize',10)
    xlabel('mean latency (ms)'),ylabel('variance latency')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),2),my_output(~isnan(my_output(:,2)),3),'Type','Pearson');
    %     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
    %     xlim([0 60])%, ylim([0 100])
    
    p_all(outInd) = p;
    r_all(outInd) = r;
    % Yeatman et al., fit a second order polynomial:
    % y  = w1* age^2 * w2*age + w3
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),2),my_output(~isnan(my_output(:,2)),3),1);
    x_mean = 1:1:100;
    y_fit = P(1)*x_mean + P(2);
    plot(x_mean,y_fit,'r')

    xlim([0 100])
    ylim([0 max(my_output(:,3))])
end

% FDR correction
m = length(p_all);
[p_sort,p_ind] = sort(p_all(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end
% figure,hold on,plot(thisVal),plot(p_sort,'r.'),title('Significant p-values after FDR correction')

% add significant stars indicating which subplots showed significant
% results after FDR corection
p_sig = p_all;
p_sig(p_ind) = p_sort<thisVal;
for outInd = 1:size(conn_matrix,1)
    subplot(4,4,outInd),hold on
    if p_sig(outInd)==1 % significant!
        plot(100,0,'r*')
    end
    
end

%% age vs var
% calculate correlation and p
p_all = zeros(size(conn_matrix,1),1);
r_all = zeros(size(conn_matrix,1),1);

figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),3);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        if length(out{outInd}.sub(kk).latencies) > 1
            my_output(kk,2) = mean(out{outInd}.sub(kk).latencies,'omitnan');
            my_output(kk,3) = var(out{outInd}.sub(kk).latencies,'omitnan');
        end
    end

    % age vs mean CCEP
    subplot(4,4,outInd),hold on
    plot(my_output(:,1),1000*my_output(:,3),'k.','MarkerSize',10)
    xlabel('age (years)'),ylabel('variance latency')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),3),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
    xlim([0 60])%, ylim([0 100])
    
    r_all(outInd,1) = r;
    p_all(outInd,1) = p;
    
    % Let's fit a first order polynomial:
    % y  =  w1*age + w2
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),1),1000*my_output(~isnan(my_output(:,2)),3),1);
    x_age = 1:1:60;
    y_fit = P(1)*x_age + P(2);
    plot(x_age,y_fit,'r')
end

% FDR correction
m = length(p_all);
[p_sort,p_ind] = sort(p_all(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end
% figure,hold on,plot(thisVal),plot(p_sort,'r.'),title('Significant p-values after FDR correction')

% add significant stars indicating which subplots showed significant
% results after FDR corection
p_sig = p_all;
p_sig(p_ind) = p_sort<thisVal;
for outInd = 1:size(conn_matrix,1)
    subplot(4,4,outInd),hold on
    if p_sig(outInd)==1 % significant!
        plot(60,0,'r*')
    end
    
end


%% age vs fano
% calculate correlation and p
p_all = zeros(size(conn_matrix,1),1);
r_all = zeros(size(conn_matrix,1),1);

figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    my_output = NaN(length(out{outInd}.sub),4);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        if length(out{outInd}.sub(kk).latencies) > 1
            my_output(kk,2) = mean(out{outInd}.sub(kk).latencies,'omitnan');
            my_output(kk,3) = var(out{outInd}.sub(kk).latencies,'omitnan');
            my_output(kk,4) = my_output(kk,3)/my_output(kk,2); % fano: variance/mean
        end
    end

    % age vs mean CCEP
    subplot(4,4,outInd),hold on
    plot(my_output(:,1),my_output(:,4),'k.','MarkerSize',10)
    xlabel('age'),ylabel('fano')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),4),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
%     xlim([0 60])%, ylim([0 100])
    p_all(outInd,1) = p;
    r_all(outInd,1) = r;
    
    % Yeatman et al., fit a second order polynomial:
    % y  = w1* age^2 * w2*age + w3
    [P,S] = polyfit(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),4),1);
    x_mean = 1:1:50;
    y_fit = P(1)*x_mean + P(2);
    plot(x_mean,y_fit,'r')

end

% FDR correction
m = length(p_all);
[p_sort,p_ind] = sort(p_all(:));
thisVal = NaN(size(p_sort));
for kk = 1:length(p_sort)
    thisVal(kk) = (kk/m)*0.05;
end
% figure,hold on,plot(thisVal),plot(p_sort,'r.'),title('Significant p-values after FDR correction')

% add significant stars indicating which subplots showed significant
% results after FDR corection
p_sig = p_all;
p_sig(p_ind) = p_sort<thisVal;
for outInd = 1:size(conn_matrix,1)
    subplot(4,4,outInd),hold on
    if p_sig(outInd)==1 % significant!
        plot(100,0,'r*')
    end
    
end