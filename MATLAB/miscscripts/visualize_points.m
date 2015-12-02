loadExperimentParams;
params.SAMPLE_NUM = 7;
params.MOVING_RUN = 7;
% LOAD IMG VOLUME 
filename = fullfile(params.INPUTDIR,sprintf('sample%dround%d_%s.tif',params.SAMPLE_NUM ,params.MOVING_RUN,params.REGISTERCHANNEL));
img = load3DTif(filename);

% LOAD KEYS

output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_sample%dround%d.mat',params.SAMPLE_NUM,params.MOVING_RUN));


load(output_keys_filename);

% figure(1);
% plot3(keyF_total(:,1),keyF_total(:,2),keyF_total(:,3),'o');
% hold on;
% for k_idx=1:size(keyF_total,1)
%     plot3(keyM_total(k_idx,1),keyM_total(k_idx,2),keyM_total(k_idx,3),'ro');
%     lines = [ ...
%             [keyM_total(k_idx,1);keyF_total(k_idx,1)] ... 
%             [keyM_total(k_idx,2);keyF_total(k_idx,2)] ...
%             [keyM_total(k_idx,3);keyF_total(k_idx,3)] ];
%      
%      rgb = [0 0 0];
%      if lines(1,1) > lines(2,1)
%          rgb(1) = .7;
%      end
%      if lines(1,2) > lines(2,2)
%          rgb(2) = .7;
%      end
%      if lines(1,3) > lines(2,3)
%          rgb(3) = .7;
%      end
%      plot3(lines(:,1),lines(:,2),lines(:,3),'color',rgb);   
% end
% [X,Y] = meshgrid(1:size(img,1),1:size(img,2));
% Z = zeros(size(X));
% frame =  im2uint8(img(:,:,30));
% % C = interp2(X, Y, P, X2, Y2);
% surf(X, Y, Z, frame, 'FaceColor', 'texturemap')
% view(90,90);

%2D Plotting:


experiments = ExperimentsManager(7);

figure(1);
subplot(4,4,1);
for i=1:length(experiments)
    subplot(4,4,i);
    exp_idx = experiments(i);
    
    if exp_idx==4; continue; end
    
    params.MOVING_RUN = exp_idx;
    output_keys_filename = fullfile(params.OUTPUTDIR,sprintf('globalkeys_sample%dround%d.mat',params.SAMPLE_NUM,params.MOVING_RUN));

    load(output_keys_filename);

    figure(1);
    plot(keyF_total(:,1),keyF_total(:,2),'o');
    hold on;
    for k_idx=1:size(keyF_total,1)
        plot(keyM_total(k_idx,1),keyM_total(k_idx,2),'ro');
        lines = [ ...
                [keyM_total(k_idx,1);keyF_total(k_idx,1)] ... 
                [keyM_total(k_idx,2);keyF_total(k_idx,2)] ...
                ];

         rgb = [0 0 0];
         if lines(1,1) > lines(2,1)
             rgb(1) = .7;
         end
         if lines(1,2) > lines(2,2)
             rgb(2) = .7;
         end

         if keyM_total(k_idx,3) ~= keyF_total(k_idx,3)
            disp(num2str(exp_idx));
         end
         plot(lines(:,1),lines(:,2),'color',rgb);   
    end
    title(['Experiment idx ' num2str(exp_idx) ]);
    hold off;

end