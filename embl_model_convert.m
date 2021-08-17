% Script to convert all EMBL models of the Machado collection
% from .xml format to .mat format 


% Function to get all files with a .xml suffix in all subdirectories
xml_files_path = uigetdir('/home/haris/Documents/coding/github_repos/GeomScale/dingo/embl_gems-master/models/');
mat_files_output_path = '/home/haris/Documents/coding/github_repos/GeomScale/dingo/embl_gems-master/mat_models/';
theFiles = dir(fullfile(xml_files_path,'**','*.xml.gz'));


for k = 1 : length(theFiles)

    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    
    gunzip(fullFileName)
    
    unzipBaseFileName = strrep(baseFileName,'.gz','');
    unzipFullFileName = strrep(fullFileName,'.gz','');
    
%   Read xml file and keep data needed
    model = readCbModel(unzipFullFileName);
    fldnames = fieldnames(model);
    dingo_model = struct;
    dingo_model.S = model.S;
    dingo_model.lb = model.lb;
    dingo_model.ub = model.ub;
    dingo_model.c = model.c;
    dingo_model.index_obj = find(dingo_model.c == 1);
    dingo_model.rxns = model.rxns;
    dingo_model.mets = model.mets;
 
    
    koo = strrep(unzipBaseFileName,'.xml','.mat');
    kico = strcat(mat_files_output_path, koo);
    fprintf(kico + "\n")

    save(kico, 'dingo_model');
    gzip(kico);
    
    delete(unzipFullFileName);
    delete(kico);

end



