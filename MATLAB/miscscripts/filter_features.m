loadExperimentParams;

% LOAD IMG VOLUME
filename = fullfile(params.INPUTDIR,sprintf('%sround%d_%s.tif',...
    params.SAMPLE_NAME,params.MOVING_RUN,params.REGISTERCHANNEL));
img = load3DTif(filename);

% LOAD KEYS

output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_%sround%d.mat',params.SAMPLE_NAME,params.MOVING_RUN));


load(output_keys_filename);
%% Create  plot of all the feature + correspondences
figure(1);

plot3(keyF_total(:,1),keyF_total(:,2),keyF_total(:,3),'o');
hold on;
for k_idx=1:size(keyF_total,1)
    plot3(keyM_total(k_idx,1),keyM_total(k_idx,2),keyM_total(k_idx,3),'ro');
    lines = [ ...
        [keyM_total(k_idx,1);keyF_total(k_idx,1)] ...
        [keyM_total(k_idx,2);keyF_total(k_idx,2)] ...
        [keyM_total(k_idx,3);keyF_total(k_idx,3)] ];
    
    rgb = [0 0 0];
    if lines(1,1) > lines(2,1)
        rgb(1) = .7;
    end
    if lines(1,2) > lines(2,2)
        rgb(2) = .7;
    end
    if lines(1,3) > lines(2,3)
        rgb(3) = .7;
    end
    plot3(lines(:,1),lines(:,2),lines(:,3),'color',rgb);
end
legend('Fixed', 'Moving');
title('Complete set of correspondences');
view(90,90);

%%
for i=1:size(keyF_total,1)
    distances(i) = sqrt( (keyM_total(i,1) - keyF_total(i,1))^2 + ...
        (keyM_total(i,2) - keyF_total(i,2))^2 + ...
        (keyM_total(i,3) - keyF_total(i,3))^2);
    
end

threshold = mean(distances);

%%

ctr_filter = 1;
for i=1:size(keyF_total,1)
    distance = sqrt( (keyM_total(i,1) - keyF_total(i,1))^2 + ...
        (keyM_total(i,2) - keyF_total(i,2))^2 + ...
        (keyM_total(i,3) - keyF_total(i,3))^2);
    
    if distance < threshold
        keyM_filtered(ctr_filter,:) = keyM_total(i,:);
        keyF_filtered(ctr_filter,:) = keyF_total(i,:);
        ctr_filter = ctr_filter +1;
    end
end

figure(2);

plot3(keyF_filtered(:,1),keyF_filtered(:,2),keyF_filtered(:,3),'o');
hold on;
for k_idx=1:size(keyF_filtered,1)
    plot3(keyM_filtered(k_idx,1),keyM_filtered(k_idx,2),keyM_filtered(k_idx,3),'ro');
    lines = [ ...
        [keyM_filtered(k_idx,1);keyF_filtered(k_idx,1)] ...
        [keyM_filtered(k_idx,2);keyF_filtered(k_idx,2)] ...
        [keyM_filtered(k_idx,3);keyF_filtered(k_idx,3)] ];
    
    rgb = [0 0 0];
    if lines(1,1) > lines(2,1)
        rgb(1) = .7;
    end
    if lines(1,2) > lines(2,2)
        rgb(2) = .7;
    end
    if lines(1,3) > lines(2,3)
        rgb(3) = .7;
    end
    plot3(lines(:,1),lines(:,2),lines(:,3),'color',rgb);
end
legend('Fixed', 'Moving');
title('Filtered set of correspondences');
view(90,90);
