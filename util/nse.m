function [ nse ] = NSE( obs,mod )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nse = 1 - (nansum((obs - mod).^2))./nansum((obs - ...
    nanmean(mod)).^2);
end

