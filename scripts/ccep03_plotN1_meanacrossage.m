


%%
%% ********* This script is work in progress *********
%%


%% load the n1Latencies from the derivatives

myDataPath = setLocalDataPath(1);
if exist(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'file')
    
    % if the n1Latencies_V1.mat was saved after ccep02_loadN1, load the n1Latencies structure here
    load(fullfile(myDataPath.output,'derivatives','av_ccep','n1Latencies_V1.mat'),'n1Latencies')
else
    disp('Run first ccep02_loadN1.mat')
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


%% try average per age first


cod_out = zeros(size(conn_matrix,1),2);
figure('position',[0 0 700 600])
for outInd = 1:size(conn_matrix,1)

    % initialize output: age, mean and variance in latency per subject
    nsubs = length(out{outInd}.sub);
    my_output = NaN(nsubs,3);

    % get variable per subject
    for kk = 1:nsubs
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = nanmean(out{outInd}.sub(kk).latencies);
        my_output(kk,3) = nanstd(out{outInd}.sub(kk).latencies)./sqrt(length(out{outInd}.sub(kk).latencies));
    end
    
    % average per age 
    x_vals = unique(sort(my_output(~isnan(my_output(:,2)),1))); % age for ~isnan
    y_vals = zeros(size(x_vals)); % this is where we will put the average for each age
    for kk = 1:length(x_vals) % each age
        y_vals(kk) = 1000*mean(my_output(ismember(my_output(:,1),x_vals(kk)),2),'omitnan');
    end
    
    % Test fitting a first order polynomial (leave 1 out cross validation)
    % y  =  w1*age + w2
    cross_val_linear = NaN(length(y_vals),4);
    % size latency (ms) X prediction (ms) X p1 (slope) X p2 (intercept) of left out
    sub_counter = 0;
    for kk = 1:length(y_vals)
        sub_counter = sub_counter+1;
        % leave out kk
        theseSubsTrain = ~ismember(1:length(y_vals),kk)';
        P = polyfit(x_vals(theseSubsTrain),y_vals(theseSubsTrain),1);
        cross_val_linear(sub_counter,3:4) = P;
        cross_val_linear(sub_counter,1) = y_vals(kk); % kk (left out) actual
        cross_val_linear(sub_counter,2) = P(1)*x_vals(kk)+P(2); % kk (left out) prediction
    end
    cod_out(outInd,1) = calccod(cross_val_linear(:,2),cross_val_linear(:,1),1);
    
    % Like Yeatman et al., for DTI fit a second order polynomial:
    cross_val_second = NaN(length(y_vals),5);
    % size latency (ms) X prediction (ms) X p1 (age^2) X p2 (age) X p3 (intercept) of left out
    sub_counter = 0;
    for kk = 1:length(y_vals)
        sub_counter = sub_counter+1;
        % leave out kk
        theseSubsTrain = ~ismember(1:length(y_vals),kk)';
        P = polyfit(x_vals(theseSubsTrain),y_vals(theseSubsTrain),2);
        cross_val_second(sub_counter,3:5) = P;
        cross_val_second(sub_counter,1) = y_vals(kk);
        cross_val_second(sub_counter,2) = P(1)*x_vals(kk).^2+P(2)*x_vals(kk)+P(3);
    end
    cod_out(outInd,2) = calccod(cross_val_second(:,2),cross_val_second(:,1),1);
    
    % Fit with a piecewise linear model:
    cross_val_piecewiselin = NaN(length(y_vals),5);
    % size latency (ms) X prediction (ms) X p1 (age^2) X p2 (age) X p3 (intercept) of left out
    sub_counter = 0;
    my_options = optimoptions(@lsqnonlin,'Display','off','Algorithm','trust-region-reflective');
    for kk = 1:length(y_vals)
        sub_counter = sub_counter+1;
        % leave out kk
        theseSubsTrain = ~ismember(1:length(y_vals),kk)';

%             % use the slmengine tool from here:
%             % John D'Errico (2020). SLM - Shape Language Modeling (https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling), MATLAB Central File Exchange. Retrieved July 14, 2020.
%             slm = slmengine(my_output(theseSubsTrain,1),1000*my_output(theseSubsTrain,2),'degree',1,'plot','off','knots',3,'interiorknots','free');
%             % predicted value at left out x
%             yhat = slmeval(my_output(kk,1),slm,0);
%             cross_val_piecewiselin(sub_counter,1) = 1000*my_output(kk,2);
%             cross_val_piecewiselin(sub_counter,2) = yhat;

        % use our own function:
        x = x_vals(theseSubsTrain);
        y = y_vals(theseSubsTrain);
        [pp] = lsqnonlin(@(pp) ccep_fitpiecewiselinear(pp,y,x),...
            [40 -1 0 20],[0 -Inf -Inf 10],[40 0 Inf 30],my_options);

        x_fit = x_vals(kk);
        y_fit = (pp(1) + pp(2)*min(pp(4),x_fit) + pp(3)*max(pp(4),x_fit));

        cross_val_piecewiselin(sub_counter,1) = y_vals(kk);
        cross_val_piecewiselin(sub_counter,2) = y_fit;
    end
    cod_out(outInd,3) = calccod(cross_val_piecewiselin(:,2),cross_val_piecewiselin(:,1),1);
    
    subplot(4,4,outInd),hold on

    % figure,hold on
    for kk = 1:nsubs
        % plot histogram per subject in the background
        if ~isnan(my_output(kk,2))
            distributionPlot(1000*out{outInd}.sub(kk).latencies','xValues',out{outInd}.sub(kk).age,...
                'color',[.6 .6 .6],'showMM',0,'histOpt',2)
        end
%         % plot mean+sterr per subject
%         plot([my_output(kk,1) my_output(kk,1)],[1000*(my_output(kk,2)-my_output(kk,3)) 1000*(my_output(kk,2)+my_output(kk,3))],...
%             'k','LineWidth',1)
    end
    % plot mean per subject in a dot
    plot(my_output(:,1),1000*my_output(:,2),'ko','MarkerSize',6)

    plot(x_vals,y_vals,'r.','MarkerSize',6)
    
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
%     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    
    % plot confidence intervals for linear fit
    x_age = [1:1:51];
    % get my crossval y distribution
    y_n1LatCross = cross_val_linear(:,3)*x_age + cross_val_linear(:,4);
    % y_n1LatCross = cross_val_second(:,3)*x_age.^2 + cross_val_second(:,4)*x_age + cross_val_second(:,5);
    
    % get 95% confidence intervals
    low_ci = quantile(y_n1LatCross,.025,1);
    up_ci = quantile(y_n1LatCross,.975,1);
    fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[0 .7 1],'EdgeColor',[0 .7 1])
    
    % put COD in title
    title(['COD=' int2str(cod_out(outInd,1)) ' p=' num2str(p,2)])
    
    xlim([0 60]), ylim([0 100])

%     xlabel('age (years)'),ylabel('mean dT (ms)')
    set(gca,'XTick',10:10:50,'YTick',20:20:100,'FontName','Arial','FontSize',12)
    
end

if ~exist(fullfile(myDataPath.output,'derivatives','age'),'dir')
    mkdir(fullfile(myDataPath.output,'derivatives','age'));    
end
figureName = fullfile(myDataPath.output,'derivatives','age','AgeVsLatency_N1_meanacrossage');


