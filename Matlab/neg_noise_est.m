function distr = neg_noise_est

rng('default')
% Having a look at the negative control probes of the lung cancer data

% These are the 770 negative control probes found on each chip
allnegCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/neg_ctrls.csv',1,1);
% labnr = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/track_labnr.csv',1,1);

n_probes = size(allnegCtrl,1)
probe_nr = (1:n_probes)';
% Might be a good idea not to use all of them. I'll leave out one third third.
probenr_test = randsample(n_probes,floor(n_probes/3));
size(probenr_test)
csvwrite('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/probenr_test.csv',probenr_test)


probenr_full = setdiff(probe_nr,probenr_test);
negCtrl = allnegCtrl(probenr_full,:);
n_probes = size(negCtrl,1)
n_ladies = size(negCtrl,2)
display_lads = random('unid',n_ladies,[1 2])

% These are the 47' regular probes found on each chip
% regular = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/regular_probes.csv',1,1);
% size(regular)

% k = [7 7 8]; n_probes = 770
k = [8 9 10]; % n_probes = 554
probenr_neg = probenr_full;

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
  probenr_neg(N10>k(j)) = [];
  % If the propability having a specific probe among the top 10 in
  % k/n_ladies trials is very low (< 0.05), the neg ctrl is excluded. k is
  % calculated by 1-binocdf(k,N,p) = 0.05, where p = 10/n_probes, N = n_ladies; 
  
  
end
save_to_base(1)
probenr_crosshybr = setdiff(probenr_full,probenr_neg);
size(probenr_crosshybr)
csvwrite('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/probenr_crosshybr.csv',probenr_crosshybr)

cross_hybr(allnegCtrl(probenr_crosshybr,:),0, display_lads,j);



