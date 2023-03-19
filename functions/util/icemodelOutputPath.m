function pathout = icemodelOutputPath(sitename,simmodel,userdata,varargin)
%ICEMODELOUTPUTPATH build icemodel output path

if nargin == 4 % fourth input is simyears
   pathout = fullfile(getenv('ICEMODELOUTPUTPATH'),sitename,simmodel,userdata, ...
      tocolumn(string(simyears)));
else
   pathout = fullfile(getenv('ICEMODELOUTPUTPATH'),sitename,simmodel,userdata);
end