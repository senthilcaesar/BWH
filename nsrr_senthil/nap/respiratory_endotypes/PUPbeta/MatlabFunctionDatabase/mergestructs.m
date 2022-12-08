function A = mergestructs(A,B,overwrite)
if ~exist('overwrite')
    overwrite=1;
end

switch overwrite
    case 1
        %data in B -> A, will overwrite existing data in A 
            f = fieldnames(B);
            for i = 1:length(f)
                A.(f{i}) = B.(f{i});
            end
    case 0
        %data in B -> A, will not overwrite existing data in A 
        a = fieldnames(A);
        b = fieldnames(B);
        I = find(sum(string(a)==string(b)')'==0); %find when they do not match
        for currField = b(I)'
            A = setfield(A,currField{:},getfield(B,currField{:}));
        end
end


