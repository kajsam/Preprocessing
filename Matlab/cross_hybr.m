function N10 = cross_hybr(negCtrl,regular)

reg_probe = 0 

% Taking a look at the means of probes
mean_negCtrl = mean(negCtrl,2);     % mean of each neg ctrl
meen_negCtrl = mean(mean_negCtrl)   % total mean of neg ctrl
max_mean_negCtrl = max(mean_negCtrl) % maximum of mean of each neg ctrl

mean_regular = mean(regular,2);     % mean of each regular probe
[sort_mean_regular,mu_idx] = sort(mean_regular); % sorted

mu_start  = find(sort_mean_regular > meen_negCtrl,1);

% having a look
[mean_regular(mu_idx(mu_start:mu_start+9)) mean_negCtrl(1:10)];

% Parameters for the EM-algorithm:  
maxiter = 500;  % Maximum number of iterations.
reps = 5; % Number of repetitions. For each rep, a new starting point.
reg = 1e-6; % Avoiding non-invertible matrices
probtol = 1e-6; % Stopping critetion
start = 'randSample'; % Starting points. 
covtype = 'full';
options = statset('MaxIter',maxiter); 

K_min = 1; K_max = 5 % Which number of components to check

n_ladies = size(negCtrl,2)
n_probes = size(negCtrl,1)
dis = zeros(2,n_ladies);

max10_neg_probes_idx = zeros(10*n_ladies);
med_neg_probes_idx = zeros(10*n_ladies);
min_neg_probes_idx = zeros(10*n_ladies);


max5_neg_probes_idx = zeros(5*n_ladies);
max1_neg_probes_idx = zeros(1*n_ladies);
t = 0;
m = 10;
for j = 1: n_ladies
  
  x = negCtrl(:,j);  % neg ctrl values for lady nr j
  [~,xidx] = sort(x); % sort them

  % Which neg ctrl probes have high values?
  max10_neg_probes_idx((j-1)*m+1:m*j) = xidx(end-(m-1):end);
  % median values?
  med_neg_probes_idx((j-1)*m+1:m*j) = xidx(floor(n_probes/2)-4:floor(n_probes/2)+5);
  % low values?
  min_neg_probes_idx((j-1)*m+1:m*j) = xidx(1:m);
  
  % Which neg ctrl probes have 5 high values?
  max5_neg_probes_idx((j-1)*5+1:5*j) = xidx(end-(5-1):end);
  % highest value?
  max1_neg_probes_idx(j) = xidx(end);
  
  if mod(j,100)==0
      t = t+1
    
    % Fitting a GMM to neg ctrls
    fit_distr = cell(1,K_max-K_min+1);
    BIC = zeros(1,K_max-K_min+1);
    AIC = zeros(1,K_max-K_min+1);
    for k = K_min: K_max
      fit_distr{k} = fitgmdist(x,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                              'Start',start,... % 'ProbabilityTolerance',probtol,
                               'CovarianceType',covtype);
      BIC(k) = fit_distr{k}.BIC;
      AIC(k) = fit_distr{k}.AIC;
    end
  
    k_BIC = BIC==min(BIC(K_min:K_max));
    k_AIC = AIC==min(AIC(K_min:K_max));
    distr = fit_distr{k_BIC};
  
    [~,idxC] = sort(distr.PComponents);
    [~,idxM] = sort(distr.mu,'descend');
    sigm = squeeze(distr.Sigma)';
  
    kBIC = find(k_BIC);
    kAIC = find(k_AIC);
    dis(:,j) = [kBIC min(distr.PComponents)]';
  
%     % If the smallest component also has the highest mean
%     if idxC(1)==idxM(1)
%       [distr.PComponents(idxC);
%        distr.mu(idxC)';
%        sigm(idxC)]
%     else % Tell me if not
%       [distr.PComponents;
%        distr.mu';
%        sigm]
%       j
%     end
      
    % Display normalised histogram and fitted GMM
    pdf_plot = pdf(distr,(min(x):0.1:max(x))');
    figure(6), histogram(x,'Normalization','pdf'), hold on
    plot((min(x):0.1:max(x))',pdf_plot)
    %subplot(1,3,2), hold on, plot(BIC), plot(AIC,'r')
    xlabel([kBIC kAIC]')
    
    figure(6+t), histogram(x,'Normalization','pdf'), hold on
    plot((min(x):0.1:max(x))',pdf_plot)
    %subplot(1,3,2), hold on, plot(BIC), plot(AIC,'r')
    xlabel([kBIC kAIC]')
    
  end
  
  if reg_probe
    %%%%%  Same procedure for n_probes of the regular probes 
    dis_reg = zeros(2,n_ladies);
    start_idx = 45001;
    if mod(j,10)==0
      x = regular(mu_idx(start_idx:start_idx+n_probes),j);
  
      fit_distr = cell(1,K_max-K_min+1);
      BIC = zeros(1,K_max-K_min+1);
      for k = K_min: K_max
        fit_distr{k} = fitgmdist(x,k,'Regularize',reg,'Options',options,'Replicates',reps, ...
                                 'Start',start,... % 'ProbabilityTolerance',probtol,
                                 'CovarianceType',covtype);
        BIC(k) = fit_distr{k}.BIC;
      end
  
      k_new = BIC==min(BIC(K_min:K_max));
      distr = fit_distr{k_new};
  
      
      [~,idxM] = sort(distr.mu,'descend');
      sigm = squeeze(distr.Sigma)';
  
      [distr.PComponents(idxM);
       distr.mu(idxM)';
       sigm(idxM)];
        
      k = find(k_new);
      dis_reg(:,j) = [k min(distr.PComponents)]';
  
      pdf_plot = pdf(distr,(min(x):0.1:max(x))');
      subplot(1,3,2), histogram(x,'Normalization','pdf'), hold on
      plot((min(x):0.1:max(x))',pdf_plot)
    end
  
    %%%%%  Same procedure for log of the regular probes
  
    x = regular(mu_idx(start_idx:start_idx+n_probes),j);
    logx = log(x);
  
    if mod(j,10)==0
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
  
      [~,idxM] = sort(distr.mu,'descend');
      sigm = squeeze(distr.Sigma)';
  
      [distr.PComponents(idxM);
       distr.mu(idxM)';
       sigm(idxM)];
          
      pdf_plot = pdf(distr,(min(logx):0.1:max(logx))');
      subplot(1,3,3), histogram(logx,'Normalization','pdf'), hold on
      plot((min(logx):0.1:max(logx))',pdf_plot)
      % pause
    end
  end
end

% Which neg ctrl probes have high values?
  figure(2), subplot(1,3,1), histogram(max10_neg_probes_idx,1:n_probes), xlabel('High 10')
  % median values?
  figure(2), subplot(1,3,2), histogram(med_neg_probes_idx,1:n_probes), xlabel('Median')
  % low values?
  figure(2), subplot(1,3,3), histogram(min_neg_probes_idx,1:n_probes), xlabel('Low')
  
  % Which neg ctrl probes have 10 high values?
  figure(3), subplot(1,3,1), histogram(max10_neg_probes_idx,1:n_probes), xlabel('High 10')
  [N10,~] = histcounts(max10_neg_probes_idx,1:n_probes);
  [~,idx10] = sort(N10,'descend');
  % 5 high values?
  figure(3), subplot(1,3,2), histogram(max5_neg_probes_idx,1:n_probes), xlabel('High 5')
  [N5,~] = histcounts(max5_neg_probes_idx,1:n_probes);
  [~,idx5] = sort(N5,'descend');
  % highest value?
  figure(3), subplot(1,3,3), histogram(max1_neg_probes_idx,1:n_probes), xlabel('High 1')
  [N1,~] = histcounts(max1_neg_probes_idx,1:n_probes);
  [~,idx1] = sort(N1,'descend');
  
  % If the propability having a specific probe among the top 10 in
  % k/n_ladies trials is very low (< 0.01), the neg ctrl is excluded. k is
  % calculated by 1-binocdf(k,N,p) = 0.01, where p = 10/n_probes, N = 272; 
  
  save_to_base(1)
   
  