function distr = neg_noise_est

negCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Preprocessing/neg_ctrls.csv',1,1);

size(negCtrl)

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

for j = 1: size(negCtrl,2)
  x = negCtrl(:,j);
  
  
  
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
  [~,idxM] = sort(distr.mu','descend');
  sigm = squeeze(distr.Sigma)';
  
  
  save_to_base(1)
  
  if idxC(1)==idxM(1)
    idx = idxC;
  
    [distr.PComponents(idx);
     distr.mu(idx)';
     sigm(idx)]
  else
    [distr.PComponents;
     distr.mu';
     sigm]
 j
 
 
  end
      
  
  k = find(k_new);
  dis(:,j) = [k min(distr.PComponents)]';
  
  pdf_plot = pdf(distr,(min(x):0.1:max(x))');
  
  clf(figure(1)), histogram(x,'Normalization','pdf'), hold on
  plot((min(x):0.1:max(x))',pdf_plot)

 pause
  
end

dis'

