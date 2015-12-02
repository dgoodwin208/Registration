%both finds the keys and converts them within the bounds of the tile
%by subtracting the y_min and x_min, which is necessary for transform calcs
function relevant_keys = findRelevantKeys(keys, y_min, y_max,x_min,x_max)
    relevant_keys = {}; relevant_ctr = 1;
    for idx=1:length(keys)
        if keys{idx}.y>=y_min && keys{idx}.y<=y_max && keys{idx}.x>=x_min && keys{idx}.x<=x_max
           relevant_keys{relevant_ctr} = keys{idx};
           relevant_keys{relevant_ctr}.y =relevant_keys{relevant_ctr}.y - y_min;
           relevant_keys{relevant_ctr}.x =relevant_keys{relevant_ctr}.x - x_min;
           relevant_ctr = relevant_ctr+1;
        end
    end
end
