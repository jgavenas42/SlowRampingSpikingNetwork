function div = kldivergence(A,ref)

%% KL Divergence between two distributions, A and a reference

% Inputs:
% A: a given distribution
% ref: a reference distribution.

% Outputs
% div: the K-L divergence between the two distributions.


% Check that the distribution sizes match
if any(size(ref) ~= size(A))
    warning("Distribution sizes don't match, using uniform as reference")
    ref = ones(size(A))/length(A);
end

% Check that the distributions add up to 1, and divide so they do if not
% if sum(A,'all') ~= 1 || sum(ref,'all' ~=1)
%     warning("Normalizing distributions")
%     A = A ./ sum(A,'all');
%     ref = ref ./ sum(ref,'all');
% end
    

div = nansum(A .* log(A./ref),'all');

end


