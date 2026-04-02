function [data1, data2] = buildOutputPayload(opts, surface_state, ...
      subsurface_state, thf_diag)
   %BUILDOUTPUTPAYLOAD Return SAVEOUTPUT cell payloads in OPTS order.
   %
   %  [data1, data2] = icemodel.buildOutputPayload(opts, surface_state, ...
   %     subsurface_state)
   %  [data1, data2] = icemodel.buildOutputPayload(opts, surface_state, ...
   %     subsurface_state, thf_diag)
   %
   % This helper is the canonical bridge between the finalized output-profile
   % contract in OPTS.VARS1 / OPTS.VARS2 and the raw timestep state assembled
   % inside the core model loops.

   if nargin < 4
      thf_diag = struct([]);
   end

   surface_output = surface_state;
   surface_output = mergeStruct(surface_output, normalizeThfDiag(thf_diag));

   data1 = selectOutputFields(surface_output, opts.vars1, 'vars1');
   data2 = selectOutputFields(subsurface_state, opts.vars2, 'vars2');
end

function output = normalizeThfDiag(diag)
   %NORMALIZETHFDIAG Map scheme-specific THF diagnostics to stable scalars.

   output = struct( ...
      'thf_es_sfc', NaN, ...
      'thf_stability_factor', NaN, ...
      'thf_z0m', NaN, ...
      'thf_z0h', NaN, ...
      'thf_z0q', NaN, ...
      'thf_u_star', NaN, ...
      'thf_L', NaN, ...
      'thf_Re', NaN, ...
      'thf_numiter', NaN);

   if isempty(diag)
      return
   end

   output.thf_es_sfc = getDiagField(diag, 'es_sfc', output.thf_es_sfc);
   output.thf_stability_factor = getDiagField(diag, 'stability_factor', ...
      output.thf_stability_factor);
   output.thf_z0m = getDiagField(diag, 'z0m', output.thf_z0m);
   output.thf_z0h = getDiagField(diag, 'z0h', output.thf_z0h);
   output.thf_z0q = getDiagField(diag, 'z0q', output.thf_z0q);
   output.thf_u_star = getDiagField(diag, 'u_star', output.thf_u_star);
   output.thf_L = getDiagField(diag, 'L', output.thf_L);
   output.thf_Re = getDiagField(diag, 'Re', output.thf_Re);
   output.thf_numiter = getDiagField(diag, 'n_iterations', ...
      output.thf_numiter);
end

function value = getDiagField(diag, field_name, default_value)
   %GETDIAGFIELD Read a field from the optional THF diagnostic struct.

   value = default_value;
   if isstruct(diag) && isfield(diag, field_name)
      value = diag.(field_name);
   end
end

function output = mergeStruct(output, extra)
   %MERGESTRUCT Copy fields from EXTRA onto OUTPUT.

   fields = fieldnames(extra);
   for n = 1:numel(fields)
      output.(fields{n}) = extra.(fields{n});
   end
end

function data = selectOutputFields(state, names, which_names)
   %SELECTOUTPUTFIELDS Return a cell array matching the requested field order.

   data = cell(1, numel(names));
   for n = 1:numel(names)
      field_name = names{n};
      if ~isfield(state, field_name)
         error('icemodel:buildOutputPayload:unknownOutputField', ...
            'Requested %s field "%s" was not assembled in the save payload.', ...
            which_names, field_name);
      end
      data{n} = state.(field_name);
   end
end
