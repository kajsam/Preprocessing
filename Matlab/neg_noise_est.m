function distr = neg_noise_est

% Having a look at the negative control probes of the lung cancer data

est = 1 % Estimate GMMs
fig = 1 % Display figures
% These are the 770 negative control probes found on each chip
negCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/neg_ctrls.csv',1,1);
size(negCtrl)

% These are the 47' regular probes found on each chip
regular = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/prob_reg.csv',1,1);
size(regular)

 

N10 = cross_hybr(negCtrl,regular);

% If the propability having a specific probe among the top 10 in
% k/n_ladies trials is very low (< 0.05), the neg ctrl is excluded. k is
% calculated by 1-binocdf(k,N,p) = 0.05, where p = 10/n_probes, N = n_ladies; 

% repeat excluded the cross-hybridization probes

negCtrl(N10>7,:) = [];
size(negCtrl)

cross_hybr(negCtrl,regular);