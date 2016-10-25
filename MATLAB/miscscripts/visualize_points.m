loadExperimentParams;

% LOAD IMG VOLUME 
filename = fullfile(params.INPUTDIR,sprintf('%sround%d_%s.tif',...
    params.SAMPLE_NAME,params.MOVING_RUN,params.REGISTERCHANNEL));
img = load3DTif(filename);

% LOAD KEYS

output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_%sround%d.mat',params.SAMPLE_NAME,params.MOVING_RUN));


load(output_keys_filename);
%% Create  plot of all the feature + correspondences
figure;

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
title(sprintf('%i correspondences to calculate TPS warp',size(keyF_total,1)))
view(45,45);

%% 2D Plotting:


% experiments = ExperimentsManager(1);

% figure(1);
% subplot(params.ROWS_DESC,params.COLS_DESC,1);
% for i=1:2 %length(experiments)
%     subplot(params.ROWS_DESC,params.COLS_DESC,i);
%     exp_idx = experiments(i);
%     
% %     if exp_idx==4; continue; end
%     
%     params.MOVING_RUN = exp_idx;
%     output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_sample%dround%d.mat',params.SAMPLE_NUM,params.MOVING_RUN));
% 
%     load(output_keys_filename);
% 
%     figure(1);
%     plot(keyF_total(:,1),keyF_total(:,2),'o');
%     hold on;
%     for k_idx=1:size(keyF_total,1)
%         plot(keyM_total(k_idx,1),keyM_total(k_idx,2),'ro');
%         lines = [ ...
%                 [keyM_total(k_idx,1);keyF_total(k_idx,1)] ... 
%                 [keyM_total(k_idx,2);keyF_total(k_idx,2)] ...
%                 ];
% 
%          rgb = [0 0 0];
%          if lines(1,1) > lines(2,1)
%              rgb(1) = .7;
%          end
%          if lines(1,2) > lines(2,2)
%              rgb(2) = .7;
%          end
% 
%          if keyM_total(k_idx,3) ~= keyF_total(k_idx,3)
%             disp(num2str(exp_idx));
%          end
%          plot(lines(:,1),lines(:,2),'color',rgb);   
%     end
%     title(['Experiment idx ' num2str(exp_idx) ]);
%     hold off;
% 
% end