


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
linear_avparams = zeros(size(conn_matrix,1),2);
piecewise_avparams = zeros(size(conn_matrix,1),6); % 4 for params, and 2 for signrank across decrease

figure('position',[0 0 700 600])
for outInd = 1:size(conn_matrix,1)
    
    % initialize output: age, mean and variance in latency per subject
    nsubs = length(out{outInd}.sub);
    my_output = NaN(nsubs,4);
    
    % get variable per subject
    for kk = 1:nsubs
        my_output(kk,1) = out{outInd}.sub(kk).age;
        my_output(kk,2) = nanmean(out{outInd}.sub(kk).latencies);
        this_nr_latencies = length(out{outInd}.sub(kk).latencies(~isnan(out{outInd}.sub(kk).latencies)));
        my_output(kk,3) = nanstd(out{outInd}.sub(kk).latencies)./sqrt(this_nr_latencies);
        my_output(kk,4) = this_nr_latencies;
        clear this_nr_latencies
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
    linear_avparams(outInd,:) = mean(cross_val_linear(:,3:4));
    
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
    cross_val_piecewiselin = NaN(length(y_vals),6);
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
            [60 -1 1 20],[20 -Inf -Inf 10],[40 0 Inf 30],my_options);

        x_fit = x_vals(kk);
        y_fit = (pp(1) + pp(2)*min(pp(4),x_fit) + pp(3)*max(x_fit-pp(4),0));
        % --> intercept = pp(1)
        % --> tipping point = pp(4)
        % --> slope before tipping point = pp(2)
        % --> slope after tipping point = pp(3)
        
        cross_val_piecewiselin(sub_counter,1) = y_vals(kk);
        cross_val_piecewiselin(sub_counter,2) = y_fit;
        cross_val_piecewiselin(sub_counter,3:6) = pp;
    end
    cod_out(outInd,3) = calccod(cross_val_piecewiselin(:,2),cross_val_piecewiselin(:,1),1);
    cod_out(outInd,4) = length(y_vals); % number of subjects
    piecewise_avparams(outInd,1:4) = mean(cross_val_piecewiselin(:,3:6));
    if quantile(cross_val_piecewiselin(:,4),.925)<0 % 85% confidence interval entirely <0
        piecewise_avparams(outInd,5) = 1;
    else
        piecewise_avparams(outInd,5) = 0;
    end
    % 
    if quantile(cross_val_piecewiselin(:,5),.925)<0 % 85% confidence interval entirely <0
        piecewise_avparams(outInd,6) = 1;
    else
        piecewise_avparams(outInd,6) = 0;
    end
        
    
    subplot(4,4,outInd),hold on
    
        % add this if you want to see every single subject and the effect of
        % averaging within an age group
        for kk = 1:nsubs
            % plot histogram per subject in the background
            if ~isnan(my_output(kk,2))
                distributionPlot(1000*out{outInd}.sub(kk).latencies','xValues',out{outInd}.sub(kk).age,...
                    'color',[.8 .8 .8],'showMM',0,'histOpt',2)
            end
    %         % plot mean+sterr per subject
    %         plot([my_output(kk,1) my_output(kk,1)],[1000*(my_output(kk,2)-my_output(kk,3)) 1000*(my_output(kk,2)+my_output(kk,3))],...
    %             'k','LineWidth',1)
        end
%         % plot mean per subject in a dot
%         plot(my_output(:,1),1000*my_output(:,2),'ko','MarkerSize',6)
    
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
    %     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    
    if cod_out(outInd,1)>cod_out(outInd,3)
        % PLOT LINEAR FIT WITH CI
        x_age = [1:1:51];
        % get my crossval y distribution
        y_n1LatCross = cross_val_linear(:,3)*x_age + cross_val_linear(:,4);
        % PLOT SECOND ORDER POLYNOMIAL FIT WITH CI
        % y_n1LatCross = cross_val_second(:,3)*x_age.^2 + cross_val_second(:,4)*x_age + cross_val_second(:,5);
        % get 95% confidence intervals
        low_ci = quantile(y_n1LatCross,.025,1);
        up_ci = quantile(y_n1LatCross,.975,1);
        if cod_out(outInd,4)<20 % less than 20 subjects
            fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[.5 .7 1],'EdgeColor',[.5 .7 1])
        else
            fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[0 .2 1],'EdgeColor',[0 .2 1])
        end
        % put COD in title
        title(['COD=' int2str(cod_out(outInd,1))])
        
    else
        % PLOT PIECEWISE LINEAR FIT (from our own function)
        x_age = [1:1:51];
        y_n1LatCross = NaN(size(cross_val_piecewiselin,1),size(x_age,2));
        for rr = 1:size(cross_val_piecewiselin,1)% run across all crossvalidations to get confidence intervals
            y_n1LatCross(rr,:) = (cross_val_piecewiselin(rr,3) + cross_val_piecewiselin(rr,4)*min(cross_val_piecewiselin(rr,6),x_age) + ...
                cross_val_piecewiselin(rr,5)*max(x_age-cross_val_piecewiselin(rr,6),0));
        end
        % get 95% confidence intervals
        low_ci = quantile(y_n1LatCross,.025,1);
        up_ci = quantile(y_n1LatCross,.975,1);
        if cod_out(outInd,4)<20 % less than 20 subjects
            fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[1 .7 .5],'EdgeColor',[1 .7 .5])
        else
            % also plot the CI of the knee in the background
            low_knee_ci = quantile(cross_val_piecewiselin(:,6),.025);
            up_knee_ci = quantile(cross_val_piecewiselin(:,6),.975);
%             fill([low_knee_ci up_knee_ci up_knee_ci low_knee_ci],[80 80 0 0],[.8 .8 .8],'EdgeColor',[.8 .8 .8])
            fill([low_knee_ci up_knee_ci up_knee_ci low_knee_ci],[5 5 0 0],[1 0 0],'EdgeColor',[1 0 0])
            
            % plot the fit
            fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[1 .2 0],'EdgeColor',[1 .2 0])
        end
        % put COD in title
        title(['COD=' int2str(cod_out(outInd,3))])
        
    end
    
    plot(x_vals,y_vals,'k.','MarkerSize',6)
    
    xlim([0 60]), ylim([0 80])
    
    %     xlabel('age (years)'),ylabel('mean dT (ms)')
    set(gca,'XTick',10:10:50,'YTick',20:20:100,'FontName','Arial','FontSize',12)
    
end

if ~exist(fullfile(myDataPath.output,'derivatives','age'),'dir')
    mkdir(fullfile(myDataPath.output,'derivatives','age'));
end
figureName = fullfile(myDataPath.output,'derivatives','age','AgeVsLatency_N1_meanacrossage');

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',figureName)
% print('-depsc','-r300',figureName)

%%
% linear better:
cod_out(cod_out<0) = 0;
lin_better = (cod_out(:,1)-cod_out(:,3))>0;
% difference in R2 for linear better
cod_out(lin_better,1) - cod_out(lin_better,3)

% mean params for linear better
mean(linear_avparams(lin_better,:))

% mean params for piecewise linear better
mean(piecewise_avparams(~lin_better & cod_out(:,4)>20,1:4))

% does the first piecewise linear decrease differ from zero?
piecewise_avparams(~lin_better & cod_out(:,4)>20,5)

% does the second piecewise linear decrease differ from zero?
piecewise_avparams(~lin_better & cod_out(:,4)>20,6)

%% mean per age instead of first mean latency per subject


cod_out = zeros(size(conn_matrix,1),2);
figure('position',[0 0 700 600])
for outInd = 1:size(conn_matrix,1)
    
    % initialize output: age, mean and variance in latency per subject
    min_age = min([out{outInd}.sub(:).age]);
    max_age = max([out{outInd}.sub(:).age]);
    my_output = NaN(max_age-min_age,4);
    
    % get latencies per age 
    counter = 1;
    for kk = min_age:max_age
        subInd = find([out{outInd}.sub(:).age] == kk);
        
        my_output(counter,1) = kk; % age
        this_age_latencies = horzcat(out{outInd}.sub(subInd).latencies);
        my_output(counter,2) = nanmean(this_age_latencies);
        this_nr_latencies = length(this_age_latencies(~isnan(this_age_latencies)));
        my_output(counter,3) = nanstd(this_age_latencies)./sqrt(this_nr_latencies);
        my_output(counter,4) = this_nr_latencies;
        clear this_nr_latencies this_age_latencies
        counter = counter+1;
    end
    
    % remove ages with no subjects (NaN)
    x_vals = my_output(~isnan(my_output(:,2)),1); % age for ~isnan
    y_vals = 1000* my_output(~isnan(my_output(:,2)),2);
    
    %     % average per 2 years of age
    %     x_vals_temp = 2*floor(.5*x_vals);
    %     x_vals_new = unique(x_vals_temp);
    %     y_vals_new = zeros(size(x_vals_new));
    %     for kk = 1:length(x_vals_new)
    %         y_vals_new(kk) = mean(y_vals(ismember(x_vals_temp,x_vals_new(kk))));
    %     end
    %     x_vals = x_vals_new; y_vals = y_vals_new;
    
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
    cross_val_piecewiselin = NaN(length(y_vals),6);
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
            [60 -1 1 20],[20 -Inf -Inf 10],[40 0 Inf 30],my_options);
        
        x_fit = x_vals(kk);
        y_fit = (pp(1) + pp(2)*min(pp(4),x_fit) + pp(3)*max(x_fit-pp(4),0));
        % --> intercept = pp(1)
        % --> tipping point = pp(4)
        % --> slope before tipping point = pp(2)
        % --> slope after tipping point = pp(3)
        
        cross_val_piecewiselin(sub_counter,1) = y_vals(kk);
        cross_val_piecewiselin(sub_counter,2) = y_fit;
        cross_val_piecewiselin(sub_counter,3:6) = pp;
    end
    cod_out(outInd,3) = calccod(cross_val_piecewiselin(:,2),cross_val_piecewiselin(:,1),1);
    cod_out(outInd,4) = length(y_vals);
    
    subplot(4,4,outInd),hold on
    
    %     % add this if you want to see every single subject and the effect of
    %     % averaging within an age group
    %     for kk = 1:nsubs
    %         % plot histogram per subject in the background
    %         if ~isnan(my_output(kk,2))
    %             distributionPlot(1000*out{outInd}.sub(kk).latencies','xValues',out{outInd}.sub(kk).age,...
    %                 'color',[.6 .6 .6],'showMM',0,'histOpt',2)
    %         end
    % %         % plot mean+sterr per subject
    % %         plot([my_output(kk,1) my_output(kk,1)],[1000*(my_output(kk,2)-my_output(kk,3)) 1000*(my_output(kk,2)+my_output(kk,3))],...
    % %             'k','LineWidth',1)
    %     end
    %     % plot mean per subject in a dot
    %     plot(my_output(:,1),1000*my_output(:,2),'ko','MarkerSize',6)
    
    [r,p] = corr(my_output(~isnan(my_output(:,2)),1),my_output(~isnan(my_output(:,2)),2),'Type','Pearson');
    %     title([out(outInd).name ' to ' out(outInd).name   ', r=' num2str(r,3) ' p=' num2str(p,3)])
    
    if cod_out(outInd,1)>cod_out(outInd,3)
        % PLOT LINEAR FIT WITH CI
        x_age = [1:1:51];
        % get my crossval y distribution
        y_n1LatCross = cross_val_linear(:,3)*x_age + cross_val_linear(:,4);
        % PLOT SECOND ORDER POLYNOMIAL FIT WITH CI
        % y_n1LatCross = cross_val_second(:,3)*x_age.^2 + cross_val_second(:,4)*x_age + cross_val_second(:,5);
        % get 95% confidence intervals
        low_ci = quantile(y_n1LatCross,.025,1);
        up_ci = quantile(y_n1LatCross,.975,1);
        fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[0 .7 1],'EdgeColor',[0 .7 1])
        % put COD in title
        title(['COD=' int2str(cod_out(outInd,1))])
        
    else
        % PLOT PIECEWISE LINEAR FIT (from our own function)
        x_age = [1:1:51];
        y_n1LatCross = NaN(size(cross_val_piecewiselin,1),size(x_age,2));
        for rr = 1:size(cross_val_piecewiselin,1)% run across all crossvalidations to get confidence intervals
            y_n1LatCross(rr,:) = (cross_val_piecewiselin(rr,3) + cross_val_piecewiselin(rr,4)*min(cross_val_piecewiselin(rr,6),x_age) + ...
                cross_val_piecewiselin(rr,5)*max(x_age-cross_val_piecewiselin(rr,6),0));
        end
        % get 95% confidence intervals
        low_ci = quantile(y_n1LatCross,.025,1);
        up_ci = quantile(y_n1LatCross,.975,1);
        fill([x_age x_age(end:-1:1)],[low_ci up_ci(end:-1:1)],[1 .7 0],'EdgeColor',[1 .7 0])
        
        % put COD in title
        title(['COD=' int2str(cod_out(outInd,3))])
        
    end
    
    plot(x_vals,y_vals,'k.','MarkerSize',6)
    
    xlim([0 60]), ylim([0 80])
    
    %     xlabel('age (years)'),ylabel('mean dT (ms)')
    set(gca,'XTick',10:10:50,'YTick',20:20:100,'FontName','Arial','FontSize',12)
    
end

if ~exist(fullfile(myDataPath.output,'derivatives','age'),'dir')
    mkdir(fullfile(myDataPath.output,'derivatives','age'));
end
figureName = fullfile(myDataPath.output,'derivatives','age','AgeVsLatency_N1_meanacrossage');
