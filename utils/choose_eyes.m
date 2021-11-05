function selected = choose_eyes(subjects,eye_option)
% chooseEyes - Select eyes to be included
%
% eyes = chooseEyes(subjects)
%
% Input arguments:
%   subjects: subjects table
%   eye_option: string 'all' (include all eyse), 'one' (choose only one eye 
%   per subject. When a subject only has one eye (due to exclusion) that 
%   eye is selected, otherwise an eye is chosen randomly.
%
% Output arguments:
%   selected: table with selected eyes
%
%  2021, Mondragon Unibertsitatea, Biomedical Engineering Department


if strcmp(eye_option,'all')
    j=1;
    for i=1:size(subjects,1)
        if subjects.has_OD(i)
            subject(j) = subjects.subject(i);
            eye(j) = {'OD'};
            age(j) = subjects.age(i);
            sex(j) = subjects.sex(i);
            j = j + 1;
        end
        
        if subjects.has_OS(i)
            subject(j) = subjects.subject(i);
            eye(j) = {'OS'};
            age(j) = subjects.age(i);
            sex(j) = subjects.sex(i);
            j = j + 1;
        end
    end

    selected = table(subject',eye',age',sex','VariableNames',{'subject',...
        'eye','age','sex'});
    
elseif strcmp(eye_option,'one')
    
    eye = cell(size(subjects,1),1);

    hasOD = subjects.has_OD;
    hasOS = subjects.has_OS;
    
    % Subjects with only OS
    eye(hasOS & ~hasOD) = {'OS'};

    % Subjects with only OD
    eye(~hasOS & hasOD) = {'OD'};

    % Subjects with both eyes (choose randomly)
    rng(0); % seed
    hasBoth = find(hasOS & hasOD);
    selectedEye = randi(2,1,length(hasBoth));
    for n=1:length(hasBoth)
        if selectedEye(n) == 1
            eye(hasBoth(n)) = {'OS'};
        elseif selectedEye(n) == 2
            eye(hasBoth(n)) = {'OD'};
        end
    end

    selected = table(subjects.subject,eye,subjects.age,subjects.sex,...
        'VariableNames',{'subject','eye','age','sex'});
else
    error('Not supported');
end




