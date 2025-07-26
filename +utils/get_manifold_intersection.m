function [idx_u, idx_s] = get_manifold_intersection(Wu,Ws)
diffs_Mat = zeros(size(Wu,1),size(Ws,1));

for i = 1:size(Wu,1)
    for j = 1:size(Ws,1)

        diffs_Mat(i,j) = sqrt(0.1*(Wu(i,2)-Ws(j,2))^2 +...
            1*(Wu(i,4)-Ws(j,4))^2 + 1*(Wu(i,3)-Ws(j,3))^2);

    end
end

[~, linearIdx] = min(diffs_Mat(:));
[idx_u, idx_s] = ind2sub(size(diffs_Mat), linearIdx);
end