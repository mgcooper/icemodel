function [heights, notes] = readKtransectHeights(filename, kwargs)
   %READKTRANSECTHEIGHTS Parse the K-transect sensor-height workbook.
   %
   %  [heights, notes] = ...
   %     icemodel.forcing.helpers.readKtransectHeights(filename)
   %  [heights, notes] = ...
   %     icemodel.forcing.helpers.readKtransectHeights(filename, ...
   %     station="AWS9")
   %
   % Role
   %  Source-specific parser for the PANGAEA.947483 series attachment
   %  "metadata_GRL_AWS_56910_2003_2021.xlsx" ("Metadata GRL AWS 2003-2021
   %  with info on sensor heights"). The "sensor heights" sheet is organized
   %  in station blocks: a numeric station id (5/6/9/10) opens each block,
   %  followed by per-year rows labeled arrival/depart (a combined
   %  "arrival/leave" for 2003) carrying the wind-sensor and
   %  temperature-sensor heights and the AWS generation in effect.
   %
   % Values policy
   %  Heights are recorded exactly as published (the sheet does not label
   %  units; the magnitudes are consistent with centimeters). Confirming the
   %  unit against Smeets et al. (2022) is part of the crosswalk/heights
   %  evidence pass; no conversion happens here.
   %
   % Name-value
   %  station : string. Optional station filter ("AWS9"); default returns
   %     every station's rows.
   %
   % Returns
   %  heights : struct array with fields station, year, event
   %     ("arrival"/"depart"/"arrival/leave"), height_u, height_T, aws_type,
   %     comment — one row per workbook height entry, in sheet order.
   %  notes : string column of the workbook's per-type height-record
   %     explanations, deduplicated in first-seen order.
   %
   % See also: icemodel.forcing.buildKtransectData,
   %  icemodel.verification.setup.fetchKtransect

   arguments
      filename (1, 1) string
      kwargs.station (1, 1) string = ""
   end

   raw = readcell(filename, 'Sheet', 'sensor heights');
   n_rows = size(raw, 1);

   % Classify rows up front so both outputs are preallocated, not grown: a
   % height entry carries an event label with a numeric height_u, while the
   % interleaved per-type explanation strings document how the acoustic
   % record relates to the boom.
   is_note = false(n_rows, 1);
   is_entry = false(n_rows, 1);
   for r = 1:n_rows
      is_note(r) = isText(raw{r, 4}) ...
         && startsWith(string(raw{r, 4}), "For AWS type");
      is_entry(r) = ~is_note(r) && isText(raw{r, 4}) && isNumber(raw{r, 5});
   end
   notes = unique(reshape(string(raw(is_note, 4)), [], 1), 'stable');
   heights = repmat(entryRecord(), 1, nnz(is_entry));

   % Walk the sheet once, tracking the current station block and year.
   current_station = "";
   current_year = NaN;
   n = 0;
   for r = 1:n_rows
      if isNumber(raw{r, 1})
         % A numeric first column opens the next station block.
         current_station = "AWS" + string(raw{r, 1});
         current_year = NaN;
      end
      if isNumber(raw{r, 3})
         % The year posts on the first row of each arrival/depart pair.
         current_year = double(raw{r, 3});
      end
      if ~is_entry(r)
         continue
      end
      n = n + 1;
      heights(n) = entryRecord(current_station, current_year, raw(r, :));
   end
   heights = heights(1:n);

   if kwargs.station ~= ""
      heights = heights(string({heights.station}) == kwargs.station);
   end
end

function tf = isText(value)
   %ISTEXT True for a nonmissing text cell.
   tf = (ischar(value) || isstring(value)) && ~any(ismissing(string(value)));
end

function tf = isNumber(value)
   %ISNUMBER True for a nonmissing numeric cell.
   tf = isnumeric(value) && ~any(ismissing(value)) && ~isempty(value);
end

function record = entryRecord(station, year, row)
   %ENTRYRECORD Build one sensor-height record (prototype when no args).
   if nargin == 0
      record = struct('station', "", 'year', NaN, 'event', "", ...
         'height_u', NaN, 'height_T', NaN, 'aws_type', NaN, 'comment', "");
      return
   end
   record = entryRecord();
   record.station = station;
   record.year = year;
   record.event = string(row{4});
   record.height_u = double(row{5});
   if isNumber(row{6})
      record.height_T = double(row{6});
   end
   if isNumber(row{7})
      record.aws_type = double(row{7});
   end
   if size(row, 2) >= 8 && isText(row{8})
      record.comment = string(row{8});
   end
end
