function X = flip_coordinates(X,eye,eye_ref)
% flip_coordinates: flip X coordinates if eye is OS
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

if strcmp(eye,eye_ref)
    X = -X;
end