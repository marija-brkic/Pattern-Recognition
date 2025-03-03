function [features] = features(bin)
    [N, M] = size(bin);
    sums = sum(bin);
    [max_sum, max_ind] = max(sums);
    [row,col] = find(bin);
    col_min = min(col);
    feature1 = max_ind-col_min;
    edges = edge(bin(:,1:max_ind), 'Canny');
    feature2 = sum(sum(edges));

    %{
    
    figure();
    subplot(1,2,1);
    imshow(bin(:,1:max_ind));
    title('Izdvajanje obeležja');
    subplot(1,2,2);
    imshow(edges);
    title('Izdvajanje obeležja');

    %}
    
    features = [feature1, feature2];
end