function [idx_u, idx_s] = get_manifold_intersection(Wu,Ws,app)
diffs_Mat = zeros(size(Wu,1),size(Ws,1));

a = 0.1;
b = 1;
c = 1;
for i = 1:size(Wu,1)
    for j = 1:size(Ws,1)
        if strcmp(app.sectionChoice.Value,'eY') || strcmp(app.sectionChoice.Value,'mY')
            diffs_Mat(i,j) = sqrt(a*(Wu(i,2)-Ws(j,2))^2 +...
                b*(Wu(i,4)-Ws(j,4))^2 + c*(Wu(i,3)-Ws(j,3))^2);
        else
            diffs_Mat(i,j) = sqrt(a*(Wu(i,1)-Ws(j,1))^2 +...
                b*(Wu(i,4)-Ws(j,4))^2 + c*(Wu(i,3)-Ws(j,3))^2);
        end

    end
end

[~, linearIdx] = min(diffs_Mat(:));
[idx_u, idx_s] = ind2sub(size(diffs_Mat), linearIdx);
end