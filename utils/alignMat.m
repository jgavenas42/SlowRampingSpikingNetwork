function centered = alignMat(W,full);

% This function will take a matrix and output another matrix where the
% diagonal elements are aligned to the left, and the off-diagonal elements
% stretch out to the right from these. 

% Note that this function assumes the input matrix is symmetric.

% Things will be different depending on if the input matrix is of even or
% odd dimension. First, check if W is square.


matsize = size(W,1); % Get size of W, will be useful for loops

    
centered = NaN(size(W)); % Create empty matrix

for i = 1:matsize
    
    % Take the off-diagonal and put it on the ith column.
    centered(1:matsize-i+1,i) = diag(W,i-1); 
        
end
    
    
if full
    
    centered2 = NaN(size(W));
    
    for i=1:matsize
        % Do the same, except for the backwards diagonal
        centered2(i:matsize,i) = diag(W,1-i);
    end
    
    centered = nanmean(cat(3,centered,centered2),3);
end