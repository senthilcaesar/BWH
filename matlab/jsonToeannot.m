
%{
The mapping container addresses the scenario where an EDF and a set of annotations differ in length, but one is not sure how they should be aligned.
The PLACE command is essentially a wrapper around SOAP, which calls it multiple times, at all possible alignments of the staging and signal data
https://zzz.bwh.harvard.edu/luna/vignettes/soap-pops/#aligning-orphaned-stage-data
The values in cell 'v' are obtained from the output of PLACE module ( Please check luna_place.m code )
For every placement, it calculates the SOAP kappa, and then will, at the end, select the offset with the highest kappa as the most likely placement.
%}

k = {'119f9726-eb4c-5a0e-a7bb-9e15256149a1.json','1da3544e-dc5c-5795-adc3-f5068959211f.json','1fa6c401-d819-50f5-8146-a0bb9e2b2516.json',...
'25a6b2b0-4d09-561b-82c6-f09bb271d3be.json','3e842aa8-bcd9-521e-93a2-72124233fe2c.json','5bf0f969-304c-581e-949c-50c108f62846.json',...
'64959ac4-53b5-5868-a845-c7476e9fdf7b.json','67fa8e29-6f4d-530e-9422-bbc3aca86ed0.json','769df255-2284-50b3-8917-2155c759fbbd.json',...
'7ab8ff5f-a77f-567d-9882-f8bee0c3c9bf.json','7d778801-88e7-5086-ad1d-70f31a371876.json','a25b2296-343b-53f6-8792-ada2669d466e.json',...
'a30245e3-4a71-565f-9636-92e7d2e825fc.json','b5d5785d-87ee-5078-b9b9-aac6abd4d8de.json','bb474ab0-c2ce-573b-8acd-ef86b0fa26a2.json',...
'f2a69bdc-ed51-5e3f-b102-6b3f7d392be0.json'};

v = {[1,932],[1,955],[1,986],[1,1063],[1,620],[61,1007],[2,992],[1,1046],[1,929],[1,989],[1,999],[1,1008],[1,1122],[61,978],[1,1142],[1,960]};
mapping = containers.Map(k, v);

% Convert JSON file to eannot
hypno_map = containers.Map({'9','0','1','2','3','4'}, {'?', 'W', 'N1', 'N2', 'N3', 'R'});
data_dir = '/Users/sq566/Downloads/dreem-learning-evaluation-master/data/annot/scorer_5';
subs = readlines([data_dir '/sub.txt'], "EmptyLineRule","skip");

for idx=1:length(subs)
    disp(['Processing subject ' num2str(idx)])
    subname = subs{idx};
    filename = [data_dir '/' subname];
    annot = str2double(readlines(filename, "EmptyLineRule","skip"));
    annot = annot(all(~isnan(annot),2),:);
    annot(annot == -1) = 9;
    hypno_val = values(hypno_map, cellstr(num2str(annot)));

    if any(strcmp(k,subname))
        x = mapping(subname);
        hypno_val = hypno_val(x(1):x(2));
    end

    subname_edit = erase(subname,".json");
    out_file = [data_dir '/' subname_edit '_eannot.txt'];
    writecell(hypno_val, out_file)
    movefile(out_file, [data_dir '/' subname_edit '.eannot']);
    clear annot hypno_val;
end




