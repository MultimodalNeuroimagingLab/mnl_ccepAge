delta_latency = [...
    -0.29 1;
    -0.42 2;
    -0.54 1;
    -0.28 1;
    -0.88 2;
    -0.23 1;
    -0.42 2;
    -0.63 2;
    -0.59 2;
    -0.77 1;
    -0.62 1;
    -0.42 1;
    -0.60 1;
    -0.83 2];

[h,p,ci,stats] = ttest2(delta_latency(delta_latency(:,2)==1,1),delta_latency(delta_latency(:,2)==2,1));