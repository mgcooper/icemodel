function [met, metadata] = data2metCollection(Data, kwargs)
   %DATA2METCOLLECTION Convert one Data timetable or a cell collection to met.
   %
   %  met = icemodel.forcing.helpers.data2metCollection(Data)
   %  [met, metadata] = ... data2metCollection(Data)
   %  met = ... data2metCollection(Data, dt_out="15m", ...
   %     fillwithmissing=true)
   %
   % Applies the canonical data2met conversion and optional interval-support
   % resampling to each source timetable. A single timetable returns a single
   % timetable and metadata struct; a cell collection preserves the met-cell
   % and metadata-struct-array input shape. Final metadata records the met
   % variables and fill policy and exactly matches met.Properties.UserData.
   % This keeps Data-backed forcing builders on the same conversion and
   % finalization path without duplicating single-point and multi-point flow.
   %
   % See also: icemodel.forcing.data2met,
   %  icemodel.forcing.helpers.resampleMetTimestep

   arguments
      Data
      kwargs.dt_out (1, 1) string = ""
      kwargs.fillwithmissing (1, 1) logical = true
      kwargs.validate (1, 1) logical = true
   end

   % Normalize the input shape before conversion so every source timetable
   % follows the same operation sequence.
   collection_input = iscell(Data);
   if ~collection_input
      Data = {Data};
   end

   met = cell(size(Data));
   metadata = cell(size(Data));
   for k = 1:numel(Data)
      met{k} = icemodel.forcing.data2met(Data{k}, ...
         fillwithmissing=kwargs.fillwithmissing, ...
         validate=kwargs.validate);
      met{k} = icemodel.forcing.helpers.resampleMetTimestep( ...
         met{k}, kwargs.dt_out);

      % Finalize the public met metadata only after resampling so callers see
      % every conversion/resampling field and the exact variables delivered.
      metadata{k} = met{k}.Properties.UserData;
      if isempty(metadata{k})
         metadata{k} = struct();
      elseif ~isstruct(metadata{k}) || ~isscalar(metadata{k})
         error('icemodel:forcing:data2metCollection:invalidMetadata', ...
            'met.Properties.UserData must be empty or a scalar metadata struct')
      end
      metadata{k}.met_variables = string(met{k}.Properties.VariableNames);
      metadata{k}.fillwithmissing = kwargs.fillwithmissing;
      metadata{k} = icemodel.forcing.helpers.columnizeMetadata(metadata{k});
      met{k}.Properties.UserData = metadata{k};
   end

   % Restore the public scalar return contract; collection metadata follows the
   % same shape as its met cell array rather than being silently columnized.
   if ~collection_input
      met = met{1};
      metadata = metadata{1};
   else
      metadata = reshape([metadata{:}], size(metadata));
   end
end
