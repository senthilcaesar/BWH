% Function to rename duplicate labels
function newLabels = renameDuplicates(labels)
    % Initialize a map to keep track of label counts
    labelCount = containers.Map('KeyType', 'char', 'ValueType', 'double');
    
    % Preallocate newLabels array
    newLabels = cell(size(labels));
    
    % Loop through labels to rename duplicates
    for i = 1:length(labels)
        label = labels{i};
        if isKey(labelCount, label)
            % Increment the count for this label
            labelCount(label) = labelCount(label) + 1;
            % Create a new label with a suffix
            newLabels{i} = [label, '_', num2str(labelCount(label))];
        else
            % Initialize count for this label
            labelCount(label) = 1;
            % Use the original label
            newLabels{i} = label;
        end
    end
end