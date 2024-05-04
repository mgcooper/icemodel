function TF = ispathlike(STR)
   %ISPATHLIKE Check if STR is consistent with a path string on the current machine
   %
   %  TF = ISPATHLIKE(STR) returns TF true if STR contains the filesep character
   %  AND the USER environment variable
   %
   % Example:
   % TF = ispathlike('C:\Users\file.txt')
   % % returns true on Windows
   %
   % TF = ispathlike('/home/user/')
   % % returns true on Linux/Mac
   %
   % TF = ispathlike('/home/user/file.txt')
   % % returns true on Linux/Mac
   %
   % See also: isfullfile, isfullpath

   % Parse inputs
   narginchk(1, 2)
   STR = convertStringsToChars(STR);

   % Do a basic check that the input is text-like
   if not(ischarlike(STR))
      TF = false;
      return
   end

   % Does the string contain the filesep
   TF = ismember(filesep, STR) && notempty(fileparts(STR));

   if ispc
      TF = TF && notempty(regexp(fileparts(STR), '^[a-zA-Z]:\\', 'once'));
   else
      TF = TF && startsWith(fileparts(STR), '/');
   end

   % Don't use this, b/c not all paths begin with /Users on macos
   % elseif ismac
   %    TF = TF && startsWith(fileparts(STR), '/Users');
end

%% LICENSE

% BSD 3-Clause License
%
% Copyright (c) 2023, Matt Cooper (mgcooper)
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
