function [session_path] = remove_empty_subjects(session_path)
%REMOVE_EMPTY_SUBJECT 
%this function removes empty subjects from session_path structure. The info
%of what was removed is included as an additional .exclude field

% (2017/08/24). Empty subjects are removed and info copied to session_path.exclude  
ind = cellfun(@(x) length(x{1})==0,session_path.behavfilenames);
    session_path.exclude.subjname = session_path.subjname(ind) ;
    session_path.exclude.subjcode = session_path.subjcode(ind) ;
    session_path.exclude.sessionfilenames = session_path.sessionfilenames(ind) ;
    session_path.exclude.behavfilenames = session_path.behavfilenames(ind) ;
    session_path.subjname(ind) = [];
    session_path.subjcode(ind) = [];
    session_path.sessionfilenames(ind) = [];
    session_path.behavfilenames(ind) = [];
end

