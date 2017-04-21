function distr = neg_noise_est

rng('default')
% Having a look at the negative control probes of the lung cancer data
tic
est = 1 % Estimate GMMs
fig = 1 % Display figures
% These are the 770 negative control probes found on each chip
negCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/neg_ctrls.csv',1,1);
n_probes = size(negCtrl,1)
% Might be a good idea not to use all of them. I'll leave out one third third.
testCtrl = random('unid',n_probes,[1 floor(n_probes/3)]);
negCtrl(testCtrl,:) = [];
n_probes = size(negCtrl,1)
n_ladies = size(negCtrl,2)
display_lads = random('unid',n_ladies,[1 2])

% These are the 47' regular probes found on each chip
% regular = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/regular_probes.csv',1,1);
% size(regular)
toc
% k = [7 7 8]; n_probes = 770
k = [8 9 10]; % n_probes = 554
for j = 1: length(k)
  [N10, medABIC] = cross_hybr(negCtrl,0, display_lads,j);
  
  n_probes = size(negCtrl,1)
  p = 10/n_probes; N = n_ladies;
  prob = 1-binocdf(k(j)-1,N,p);
  [k(j)-1 prob]
  prob = 1-binocdf(k(j),N,p);
  [k(j) prob]
  prob = 1-binocdf(k(j)+1,N,p);
  [k(j)+1 prob]
  if sum(medABIC) == 2
    break
  end
 
  negCtrl(N10>k(j),:) = [];
  % If the propability having a specific probe among the top 10 in
  % k/n_ladies trials is very low (< 0.05), the neg ctrl is excluded. k is
  % calculated by 1-binocdf(k,N,p) = 0.05, where p = 10/n_probes, N = n_ladies; 
  
  %     pause
end

