function distr = neg_noise_est

% Having a look at the negative control probes of the lung cancer data

est = 1 % Estimate GMMs
fig = 1 % Display figures
% These are the 770 negative control probes found on each chip
negCtrl = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/neg_ctrls.csv',1,1);
n_probes = size(negCtrl,1)
n_ladies = size(negCtrl,2)

% These are the 47' regular probes found on each chip
regular = csvread('/Volumes/kam025/Documents/LungCancer/Discrete_curve_group_NR_method/Preprocessing/regular_probes.csv',1,1);
size(regular)

 

N10 = cross_hybr(negCtrl,regular);

% If the propability having a specific probe among the top 10 in
% k/n_ladies trials is very low (< 0.05), the neg ctrl is excluded. k is
% calculated by 1-binocdf(k,N,p) = 0.05, where p = 10/n_probes, N = n_ladies; 

p = 10/n_probes; N = n_ladies;
k = 10
1-binocdf(k,N,p)
k = 5
1-binocdf(k,N,p)
k = 7
1-binocdf(k,N,p)

% repeat excluded the cross-hybridization probes

negCtrl(N10>7,:) = [];
size(negCtrl)

cross_hybr(negCtrl,regular);