function [cv2, cv2s] = calculateCV2(spks)

%% Caculates the coefficient of variation 2 given a spike train (1D)


if size(spks,1)~=1
    spks = spks';
end

ISIs = diff(find(spks==1));

cv2_top = 2 * abs((ISIs(1:end-1) - ISIs(2:end)));
cv2_bottom = ISIs(1:end-1) + ISIs(2:end);

cv2s = cv2_top / cv2_bottom;

cv2 = mean(cv2s);
