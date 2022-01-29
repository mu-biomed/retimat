function X = flip_coordinates(X, eye, eye_ref)
% flip_coordinates: flip X coordinates if eye is equal to eye_ref
%
% function X = flip_coordinates(X,eye,eye_ref)
%
% Input arguments
% - X: input x coordinates (in any shape)
% - eye: eye to be tested
% - eye_ref: reference to decide when to flip (usually OS)
%
% Output arguments
% - X: flipped coordinates
%
% David Romero-Bascones 
% dromero@mondragon.edu
% 2021, Mondragon Unibertsitatea, Biomedical Engineering Department
% -------------------------------------------------------------------------

if ~any(strcmp(eye, {'OD','OS'})) || ~any(strcmp(eye_ref, {'OD','OS'}))
    error("The type of eye must be either OD or OS");
end

if strcmp(eye,eye_ref)
    X = -X;
end