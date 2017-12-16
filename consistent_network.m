% Version 01/12/2017
%
% Reference: Tokariev A, Stjerna S, Palva JM, Vanhatalo S. 
%           'Preterm birth changes networks of newborn cortical activity'
%
% For more details see Fig.S5 in the Supplementary to the paper


function [Cons_netw, Cons_netw_fr, D, Overlap] = consistent_network(Data, k, FDR)

% INPUT
% =====
% Data: array of cells [N subjects x 1]
%       each cell contains adjacency matrix (AM) for corresponding subject
%
%       Data = [{subjA};
%               {subjB};
%               {subjC}]
%
%       Note, AI - square, symmetric, weighted
%
%   k: Input density (k=0.1 (10%), k=0.2 (20%), ... of the strongest edges)
%
% FDR: if you want to apply FDR (Benjamini-Hochberg) correction set '1' 
%      if you do not want to apply FDR (take all p < 0.05) set '0'

% OUTPUT
% ======
%    Cons_netw: subset of edges that are consistent across group (binary)
% Cons_netw_fr: frequency of consistent edges
%            D: density of consistent network
%      Overlap: coefficient showing overlap of individual 'k-networks'
%               (or top 'k' strongest edges in a subject) with group level 
%               consistent network


% Number of subjects
  N_sb = size(Data, 1);
  
  if size(Data{1, 1}, 1) == size(Data{1, 1}, 2)
     M = size(Data{1, 1}, 1); % Number of nodes
  else
     disp('Matrix is not square! Check!');
  end
  
% Take ABS => all positive  
  Data = cellfun(@abs, Data, 'UniformOutput', false); 
  
% Set lower corner of AM to zeros
  for n = 1:N_sb
 
      buf = [];
      buf = Data{n, 1};
      buf(tril(ones(size(buf)), 0) == 1) = 0;
          
      Data{n, 1} = buf;

  end
  
% Compute number of edges in original networks
  N_edges = zeros(N_sb, 1);

  for n = 1:N_sb
      N_edges(n, 1) = length(nonzeros(Data{n, 1}));
  end  
  
  N_edges = mean(N_edges);
  
% Compute binary 'k-networks'
  Data_bin = cell(N_sb, 1);

  for n = 1:N_sb
      Data_bin{n, 1} = get_thresholded_values(Data{n, 1}, k);
  end
   
% Compute sum of binary arrays  
  Data_bin_sum = sum(cell2mat(reshape(Data_bin, 1, 1, N_sb)), 3);
  
% Binomial stat.
  Data_bin_stat = zeros(M);
  
  for chA = 1:M
     for chB = (chA + 1):M
         Data_bin_stat(chA, chB) = binocdf(Data_bin_sum(chA, chB), N_sb, k, 'upper');
     end
  end
    
% Consistent network = all edges with p < 0.05 (+FDR if FDR = '1')  
  if FDR == 1
     Cons_netw = make_FDR(Data_bin_stat);
  else
     Cons_netw = double(Data_bin_stat < 0.05 & Data_bin_stat > 0);
  end
  
% Density of consistent edges
  D = sum(sum(Cons_netw)) / N_edges;
    
% Frequency  
  Cons_netw_fr = (Data_bin_sum .* Cons_netw) ./ N_sb;
    
% Overlap ('k-networks' vs. Consistent(freqency) network)
  Cons_netw_fr_norm = Cons_netw_fr ./ sum(sum(Cons_netw_fr));
  
  Overlap = zeros(N_sb, 1); 
  for n = 1:N_sb
      Overlap(n, 1) = sum(sum(Data_bin{n, 1} .* Cons_netw_fr_norm));
  end


end



%Computes thresholds as % of strongest values
function B = get_thresholded_values(A, k)

   values = nonzeros(A);                    
    
   values = sort(values, 'descend');
   
   T = values(fix(length(values)*k));
   
   B = double(A >= T);
   
end
 


function B = make_FDR(A)

   p = nonzeros(A);                    
   
   p = sort(p, 'ascend');
   
   FDR = mafdr(p, 'BHFDR', 'true');
   
   FDR(FDR >= 0.05) = 0;
   
   [~, n] = max(FDR);
   
   T = p(n);
   
   B = double(A <= T & A > 0);

end




