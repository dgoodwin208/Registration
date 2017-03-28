% INPUTS:
% run_num_list is the index list of the experiment for the specified sample
function calculateDescriptorsInParallel(run_num_list)

    loadExperimentParams;

    run_num_list_size = length(run_num_list);
    desc_size = params.ROWS_DESC * params.COLS_DESC;
    run_size  = run_num_list_size * desc_size;

    disp('set up cluster')
    tic;
    cluster = parcluster('local_96workers');
    toc;

    tic;
    disp('create batch jobs')

    jobs = cell(1, run_size);

    for i = 1:run_size
        run_num    = run_num_list(ceil(i / desc_size));
        target_idx = mod(i, desc_size);
        if target_idx == 0
            target_idx = desc_size;
        end

        disp(['create batch (',num2str(i),') run_num=',num2str(run_num),', target_idx=',num2str(target_idx)])
        jobs{i} = batch(cluster,@calculateDescriptors,0,{run_num,target_idx,target_idx},'Pool',2,'CaptureDiary',true);

    end
    toc;

    tic;
    disp('waiting batch jobs...')
    for i = 1:run_size
        run_num    = run_num_list(ceil(i / desc_size));
        target_idx = mod(i, desc_size);
        if target_idx == 0
            target_idx = desc_size;
        end

        wait(jobs{i})
        diary(jobs{i},['./matlab-calcDesc-',num2str(run_num),'-',num2str(target_idx),'.log']);
    end

    disp('all batch jobs finished')
    toc;

    tic;
    disp('delete batch jobs')
    for i = 1:run_size
        delete(jobs{i})
    end
    toc;


end

