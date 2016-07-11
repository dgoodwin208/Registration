function [ output_args ] = dedupeOverlapFeatures(keys)
%DEDUPEOVERLAPFEATURES Each keypoint/descriptor in the cell array should
%have a unique {x,y,z,k} value. This function removes any of those duplicates 
%           x: 10
%           y: 90
%           z: 10
%     xyScale: 1
%      tScale: 1
%        ivec: [1x640 uint8]
%           k: 1.0620

%unique_tuples = zeros(1,4);
ctr = 1;
for k = 1:length(keys)
    tuple = [keys{k}.x,keys{k}.y,keys{k}.z, keys{k}.k];
    unique_tuples(ctr,:) = tuple;
    
    %%TODO: How to efficiently check for duplicates?
end

end

