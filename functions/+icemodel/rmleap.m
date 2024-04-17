function varargout = rmleap(simyear, Data, varargin)
   %RMLEAPINDS Remove leap-year indices from icemodel data structures
   % 
   % [ice1, ice2] = icemodel.rmleapinds(ice1, ice2)
   % [ice1, ice2] = icemodel.rmleapinds(ice1, ice2)
   
   % NOTE: For now this just replaces the inline code in prepRunoff where this
   % is called to simplify but it would be a good general purpose function
   
   if istimetable(varargin{1})
      ice1 = varargin{1};
   end
   
   if isleap(simyear)
      T = ice1.Time;
      Y = year(T(1));
      Tleap = datetime(Y, 1, 1, 0, 0, 0):hours(1):datetime(Y, 12, 31, 23, 0, 0);
      feb29 = month(Tleap)==2 & day(Tleap)==29;
      fields = fieldnames(Data);
      for n = 1:numel(fields)
         dat = Data.(fields{n});
         if size(dat, 1) == numel(Tleap)
            dat = rmleapinds(dat, Tleap);
            Data.(fields{n}) = dat;
         end
      end
      if numel(T) == numel(Tleap)
         ice1 = ice1(~feb29, :);
      end
   end
   
   switch nargout
      case 1
         varargout{1} = Data;
      case 2
         varargout{1} = Data;
         varargout{2} = ice1;
   end
end
