function distr = neg_noise_est

fig = 0
est = 0
% These are the 770 negative control probes found on each chip
negCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/neg_ctrls.csv',1,1);
size(negCtrl)

% These are the 47' regular probes found on each chip
probCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/prob_ctrls.csv',1,1);
size(probCtrl)

mean_negCtrl = mean(negCtrl,2);     % mean of each neg control probe
mean_probCtrl = mean(probCtrl,2);   % mean of each regular probe
[sort_mean_probCtrl,mu_idx] = sort(mean_probCtrl); % sorted

meen_negCtrl = mean(mean_negCtrl)   % total mean of neg control probes
max_meen_negCtrl = max(mean_negCtrl) % maximum of mean of each neg ctrl

mu_start  = find(sort_mean_probCtrl > meen_negCtrl,1);

% having a look
[mean_probCtrl(mu_idx(mu_start:mu_start+49)) mean_negCtrl(1:50)]; 

% Parameters for the EM-algorithm:  
maxiter = 500;  % Maximum number of iterations.
reps = 5; % Number of repetitions. For each rep, a new starting point.
reg = 1e-6; % Avoiding non-invertible matrices
probtol = 1e-6; % Stopping critetion
start = 'randSample'; % Starting points. 
covtype = 'full';
EMparam{1} = [maxiter reps reg probtol];
EMparam{2} = start;
EMparam{3} = covtype;
K_min = 1; K_max = 4;

maxiter = EMparam{1}(1);
reps = EMparam{1}(2);
reg  = EMparam{1}(3);
probtol = EMparam{1}(4);
start = EMparam{2};
covtype = EMparam{3};
options = statset('MaxIter',maxiter); 

dis = zeros(2,size(negCtrl,2));

max_neg_probes_idx = [];
med_neg_probes_idx = [];
min_neg_probes_idx = [];

max10_neg_probes_idx = [];
max5_neg_probes_idx = [];
max2_neg_probes_idx = [];

for j = 1: size(negCtrl,2)
  x = negCtrl(:,j);
  
  [~,xidx] = sort(x);
  
  % Which neg ctrl probes have high values?
  max_neg_probes_idx = [max_neg_probes_idx; xidx(end-4:end)];
  % median values?
  med_neg_probes_idx = [med_neg_probes_idx; xidx(770/2-2:770/2+2)];
  % low values?
  min_neg_probes_idx = [min_neg_probes_idx; xidx(1:5)];
  
  % Which neg ctrl probes have 10 high values?
  max10_neg_probes_idx = [max10_neg_probes_idx; xidx(end-9:end)];
  % 5 high values?
  max5_neg_probes_idx = [max5_neg_probes_idx; xidx(end-4:end)];
  % highest value?
  max2_neg_probes_idx = [max2_neg_probes_idx; xidx(end-1:end)];
  
  
  if mod(j,10)==0
     %  K_min = 1; K_max = 1;
  % Fitting a GMM to neg ctrls
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(x,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC;
  end
  
  
  k_new = BIC==min(BIC(K_min:K_max));
  BIC
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
  % If the smallest component also has the highest mean
  if idxC(1)==idxM(1)
    [distr.PComponents(idxC);
     distr.mu(idxC)';
     sigm(idxC)]
  else % Tell me if not
    [distr.PComponents;
     distr.mu';
     sigm]
    j
  end
      
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  
  % if fig
  % Display normalised histogram and fitted GMM
  pdf_plot = pdf(distr,(min(x):0.1:max(x))');
  clf(figure(1)), subplot(1,3,1), histogram(x,'Normalization','pdf'), hold on
  plot((min(x):0.1:max(x))',pdf_plot)
  pause
  % end
  end
  %%%%%  Same procedure for 770 of the regular probes
  % mu_start-769:mu_start
  
  
  x = probCtrl(mu_idx(45001:45000+770),j);
  logx = x;
  
  if est
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(logx,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC;
  end
  
  k_new = BIC==min(BIC(K_min:K_max));
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
  [distr.PComponents(idxC);
   distr.mu(idxC)';
   sigm(idxC)];
        
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  
  if fig
  pdf_plot = pdf(distr,(min(logx):0.1:max(logx))');
  subplot(1,3,2), histogram(logx,'Normalization','pdf'), hold on
  plot((min(logx):0.1:max(logx))',pdf_plot)
  end
  end
  
  %%%%%  Same procedure for 770 of the regular probes
  % mu_start-769:mu_start
  
  
  x = probCtrl(mu_idx(45001:45000+770),j);
  logx = log(x);
  
  if est
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(logx,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC;
  end
  
  k_new = BIC==min(BIC(K_min:K_max));
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
  [distr.PComponents(idxC);
   distr.mu(idxC)';
   sigm(idxC)];
        
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  if fig
  pdf_plot = pdf(distr,(min(logx):0.1:max(logx))');
  subplot(1,3,3), histogram(logx,'Normalization','pdf'), hold on
  plot((min(logx):0.1:max(logx))',pdf_plot)
  end
  end
  % pause
end

% Which neg ctrl probes have high values?
  figure(2), subplot(1,3,1), hist(max_neg_probes_idx,770), xlabel('High 5')
  % median values?
  figure(2), subplot(1,3,2), hist(med_neg_probes_idx,770), xlabel('Median')
  % low values?
  figure(2), subplot(1,3,3), hist(min_neg_probes_idx,770), xlabel('Low')
  
  n_cross = 20;
  % Which neg ctrl probes have 10 high values?
  figure(3), subplot(1,3,1), hist(max10_neg_probes_idx,1:770), xlabel('High 10')
  [N10,~] = histcounts(max10_neg_probes_idx,1:770);
  [~,idx10] = sort(N10,'descend');
    idx10(1:n_cross)
  % 5 high values?
  figure(3), subplot(1,3,2), hist(max5_neg_probes_idx,1:770), xlabel('High 5')
  [N5,~] = histcounts(max5_neg_probes_idx,1:770);
  [~,idx5] = sort(N5,'descend');
  idx5(1:n_cross)
  % highest value?
  figure(3), subplot(1,3,3), hist(max2_neg_probes_idx,1:770), xlabel('High 1')
  [N2,~] = histcounts(max2_neg_probes_idx,1:770);
  [~,idx2] = sort(N2,'descend');
  idx2(1:n_cross)
  

   idx15 = intersect(idx10(1:n_cross),idx5(1:n_cross));
   idx_cross = intersect(idx15,idx2(1:n_cross))
   
   p = 10/770; 
   N = 272; 
   
   idx_out = find(N10>7)
   %% repeat excluded the cross-hybridization probes

   fig = 0;
save_to_base(1)

negCtrl(idx_out,:) = [];
size(negCtrl)
mean_negCtrl = mean(negCtrl,2);     % mean of each neg control probe
mean_probCtrl = mean(probCtrl,2);   % mean of each regular probe
[sort_mean_probCtrl,mu_idx] = sort(mean_probCtrl); % sorted

meen_negCtrl = mean(mean_negCtrl)   % total mean of neg control probes
max_meen_negCtrl = max(mean_negCtrl) % maximum of mean of each neg ctrl

mu_start  = find(sort_mean_probCtrl > meen_negCtrl,1);

% having a look
[mean_probCtrl(mu_idx(mu_start:mu_start+49)) mean_negCtrl(1:50)]; 

% Parameters for the EM-algorithm:  
maxiter = 500;  % Maximum number of iterations.
reps = 5; % Number of repetitions. For each rep, a new starting point.
reg = 1e-6; % Avoiding non-invertible matrices
probtol = 1e-6; % Stopping critetion
start = 'randSample'; % Starting points. 
covtype = 'full';
EMparam{1} = [maxiter reps reg probtol];
EMparam{2} = start;
EMparam{3} = covtype;
K_min = 4; K_max = 4;

maxiter = EMparam{1}(1);
reps = EMparam{1}(2);
reg  = EMparam{1}(3);
probtol = EMparam{1}(4);
start = EMparam{2};
covtype = EMparam{3};
options = statset('MaxIter',maxiter); 

dis = zeros(2,size(negCtrl,2));

max_neg_probes_idx = [];
med_neg_probes_idx = [];
min_neg_probes_idx = [];

max10_neg_probes_idx = [];
max5_neg_probes_idx = [];
max2_neg_probes_idx = [];

for j = 1: size(negCtrl,2)
  x = negCtrl(:,j);
  
  [~,xidx] = sort(x);
  
  % Which neg ctrl probes have high values?
  max_neg_probes_idx = [max_neg_probes_idx; xidx(end-4:end)];
  % median values?
  med_neg_probes_idx = [med_neg_probes_idx; xidx(770/2-2:770/2+2)];
  % low values?
  min_neg_probes_idx = [min_neg_probes_idx; xidx(1:5)];
  
  % Which neg ctrl probes have 10 high values?
  max10_neg_probes_idx = [max10_neg_probes_idx; xidx(end-9:end)];
  % 5 high values?
  max5_neg_probes_idx = [max5_neg_probes_idx; xidx(end-4:end)];
  % highest value?
  max2_neg_probes_idx = [max2_neg_probes_idx; xidx(end-1:end)];
  
  
  if est
  % Fitting a GMM to neg ctrls
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(x,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC;
  end
  
  
  k_new = BIC==min(BIC(K_min:K_max));
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
%   % If the smallest component also has the highest mean
%   if idxC(1)==idxM(1)
%     [distr.PComponents(idxC);
%      distr.mu(idxC)';
%      sigm(idxC)]
%   else % Tell me if not
%     [distr.PComponents;
%      distr.mu';
%      sigm]
%     j
%   end
      
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  
  if fig
  % Display normalised histogram and fitted GMM
  pdf_plot = pdf(distr,(min(x):0.1:max(x))');
  clf(figure(1)), subplot(1,3,1), histogram(x,'Normalization','pdf'), hold on
  plot((min(x):0.1:max(x))',pdf_plot)
  end
  end
  %%%%%  Same procedure for 770 of the regular probes
  % mu_start-769:mu_start
  
  
  x = probCtrl(mu_idx(45001:45000+770),j);
  
  logx = x;
  
  if est
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(logx,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC;
  end
  
  k_new = BIC==min(BIC(K_min:K_max));
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
  [distr.PComponents(idxC);
   distr.mu(idxC)';
   sigm(idxC)];
        
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  
  if fig
  pdf_plot = pdf(distr,(min(logx):0.1:max(logx))');
  subplot(1,3,2), histogram(logx,'Normalization','pdf'), hold on
  plot((min(logx):0.1:max(logx))',pdf_plot)
  end
  end

  
  %%%%%  Same procedure for 770 of the regular probes
  % mu_start-769:mu_start
  
  
  x = probCtrl(mu_idx(45001:45000+770),j);
  logx = log(x);
  
  if est
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(logx,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC;
  end
  
  k_new = BIC==min(BIC(K_min:K_max));
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
  [distr.PComponents(idxC);
   distr.mu(idxC)';
   sigm(idxC)];
        
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  if fig
  pdf_plot = pdf(distr,(min(logx):0.1:max(logx))');
  subplot(1,3,3), histogram(logx,'Normalization','pdf'), hold on
  plot((min(logx):0.1:max(logx))',pdf_plot)
  pause
  end
  end
   
end

% Which neg ctrl probes have high values?
  figure(4), subplot(1,3,1), hist(max_neg_probes_idx,770), xlabel('High 5')
  % median values?
  figure(4), subplot(1,3,2), hist(med_neg_probes_idx,770), xlabel('Median')
  % low values?
  figure(4), subplot(1,3,3), hist(min_neg_probes_idx,770), xlabel('Low')
  
  n_cross = 20;
  % Which neg ctrl probes have 10 high values?
  figure(5), subplot(1,3,1), hist(max10_neg_probes_idx,1:770), xlabel('High 10')
  [N10,~] = histcounts(max10_neg_probes_idx,1:770);
  [~,idx10] = sort(N10,'descend');
    idx10(1:n_cross)
  % 5 high values?
  figure(5), subplot(1,3,2), hist(max5_neg_probes_idx,1:770), xlabel('High 5')
  [N5,~] = histcounts(max5_neg_probes_idx,1:770);
  [~,idx5] = sort(N5,'descend');
  idx5(1:n_cross)
  % highest value?
  figure(5), subplot(1,3,3), hist(max2_neg_probes_idx,1:770), xlabel('High 1')
  [N2,~] = histcounts(max2_neg_probes_idx,1:770);
  [~,idx2] = sort(N2,'descend');
  idx2(1:n_cross)
  

   idx15 = intersect(idx10(1:n_cross),idx5(1:n_cross));
   idx_cross2 = intersect(idx15,idx2(1:n_cross))
   
   idx_out = find(N10>9)
   %% repeat again excluded the cross-hybridization probes

   fig = 1
   est = 1
save_to_base(1)

negCtrl(idx_out,:) = [];
size(negCtrl)
mean_negCtrl = mean(negCtrl,2);     % mean of each neg control probe
mean_probCtrl = mean(probCtrl,2);   % mean of each regular probe
[sort_mean_probCtrl,mu_idx] = sort(mean_probCtrl); % sorted

meen_negCtrl = mean(mean_negCtrl)   % total mean of neg control probes
max_meen_negCtrl = max(mean_negCtrl) % maximum of mean of each neg ctrl

mu_start  = find(sort_mean_probCtrl > meen_negCtrl,1);

% having a look
[mean_probCtrl(mu_idx(mu_start:mu_start+49)) mean_negCtrl(1:50)]; 

% Parameters for the EM-algorithm:  
maxiter = 500;  % Maximum number of iterations.
reps = 5; % Number of repetitions. For each rep, a new starting point.
reg = 1e-6; % Avoiding non-invertible matrices
probtol = 1e-6; % Stopping critetion
start = 'randSample'; % Starting points. 
covtype = 'full';
EMparam{1} = [maxiter reps reg probtol];
EMparam{2} = start;
EMparam{3} = covtype;
K_min = 4; K_max = 4;

maxiter = EMparam{1}(1);
reps = EMparam{1}(2);
reg  = EMparam{1}(3);
probtol = EMparam{1}(4);
start = EMparam{2};
covtype = EMparam{3};
options = statset('MaxIter',maxiter); 

dis = zeros(2,size(negCtrl,2));

max_neg_probes_idx = [];
med_neg_probes_idx = [];
min_neg_probes_idx = [];

max10_neg_probes_idx = [];
max5_neg_probes_idx = [];
max2_neg_probes_idx = [];

for j = 1: size(negCtrl,2)
  x = negCtrl(:,j);
  
  [~,xidx] = sort(x);
  
  % Which neg ctrl probes have high values?
  max_neg_probes_idx = [max_neg_probes_idx; xidx(end-4:end)];
  % median values?
  med_neg_probes_idx = [med_neg_probes_idx; xidx(770/2-2:770/2+2)];
  % low values?
  min_neg_probes_idx = [min_neg_probes_idx; xidx(1:5)];
  
  % Which neg ctrl probes have 10 high values?
  max10_neg_probes_idx = [max10_neg_probes_idx; xidx(end-9:end)];
  % 5 high values?
  max5_neg_probes_idx = [max5_neg_probes_idx; xidx(end-4:end)];
  % highest value?
  max2_neg_probes_idx = [max2_neg_probes_idx; xidx(end-1:end)];
  
  if est && mod(j,10)==0
      K_min = 1; K_max = 1;
  % Fitting a GMM to neg ctrls
  fit_distr = cell(1,K_max-K_min+1);
  BIC = zeros(1,K_max-K_min+1);
  for k = K_min: K_max
    fit_distr{k} = fitgmdist(x,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                               'ProbabilityTolerance',probtol,'Start',start,...
                               'CovarianceType',covtype);
    BIC(k) = fit_distr{k}.BIC
  end
  
  
  k_new = BIC==min(BIC(K_min:K_max));
  distr = fit_distr{k_new};
  
  [~,idxC] = sort(distr.PComponents);
  [~,idxM] = sort(distr.mu,'descend');
  sigm = squeeze(distr.Sigma)';
  
  save_to_base(1)
  
%   % If the smallest component also has the highest mean
%   if idxC(1)==idxM(1)
%     [distr.PComponents(idxC);
%      distr.mu(idxC)';
%      sigm(idxC)]
%   else % Tell me if not
%     [distr.PComponents;
%      distr.mu';
%      sigm]
%     j
%   end
      
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  
  if fig
  % Display normalised histogram and fitted GMM
  pdf_plot = pdf(distr,(min(x):0.1:max(x))');
  clf(figure(1)), histogram(x,'Normalization','pdf'), hold on
  plot((min(x):0.1:max(x))',pdf_plot)
  pause
  end
  end
  
  
  
  
   
end

% Which neg ctrl probes have high values?
  figure(4), subplot(1,3,1), hist(max_neg_probes_idx,770), xlabel('High 5')
  % median values?
  figure(4), subplot(1,3,2), hist(med_neg_probes_idx,770), xlabel('Median')
  % low values?
  figure(4), subplot(1,3,3), hist(min_neg_probes_idx,770), xlabel('Low')
  
  n_cross = 20;
  % Which neg ctrl probes have 10 high values?
  figure(5), subplot(1,3,1), hist(max10_neg_probes_idx,1:770), xlabel('High 10')
  [N10,~] = histcounts(max10_neg_probes_idx,1:770);
  [~,idx10] = sort(N10,'descend');
    idx10(1:n_cross)
  % 5 high values?
  figure(5), subplot(1,3,2), hist(max5_neg_probes_idx,1:770), xlabel('High 5')
  [N5,~] = histcounts(max5_neg_probes_idx,1:770);
  [~,idx5] = sort(N5,'descend');
  idx5(1:n_cross)
  % highest value?
  figure(5), subplot(1,3,3), hist(max2_neg_probes_idx,1:770), xlabel('High 1')
  [N2,~] = histcounts(max2_neg_probes_idx,1:770);
  [~,idx2] = sort(N2,'descend');
  idx2(1:n_cross)
  

   idx15 = intersect(idx10(1:n_cross),idx5(1:n_cross));
   idx_cross2 = intersect(idx15,idx2(1:n_cross))
   
   save_to_base(1)