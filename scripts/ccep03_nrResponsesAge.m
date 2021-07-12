%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
end

%% age distribution

all_ages = zeros(length(n1Latencies),1);
for kk = 1:length(n1Latencies)
    all_ages(kk) = n1Latencies(kk).age;
end


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

%%
figure('position',[0 0 700 700])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean, variance in latency and #N1s per subject
    my_output = NaN(length(out{outInd}.sub),4);

    % get variable per subject
    for kk = 1:length(out{outInd}.sub)
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = mean(out{outInd}.sub(kk).latencies,'omitnan');
        my_output(kk,3) = nanvar(out{outInd}.sub(kk).latencies);
        my_output(kk,4) = numel(out{outInd}.sub(kk).latencies);
    end

    % age vs # CCEP
    subplot(4,4,outInd),hold on
    plot(my_output(:,1),my_output(:,4),'k.','MarkerSize',10)
    xlabel('age (years)'),ylabel('# CCEPs')
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),4),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    title(['r=' num2str(r,3) ' p=' num2str(p,3)])
    xlim([0 60])%, ylim([0 100])
    
end


% figureName = fullfile(myDataPath.output,'derivatives','age','AgeVsNumber_N1');

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)


