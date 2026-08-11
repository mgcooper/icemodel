function tf = artifactScalarIdentityMatches(existing, incoming)
   %ARTIFACTSCALARIDENTITYMATCHES Compare concrete scalar provenance facts.
   %
   %  tf = icemodel.forcing.helpers.artifactScalarIdentityMatches( ...
   %     existing, incoming)
   %
   % A record that omits a metadata field stays compatible. The function returns
   % false when both records know a field and the values differ. The checked
   % fields are family, source, product, native producer, relationship, DOI, and
   % schema. Production `method` and repaired or manifest `sample_method` are
   % the one documented alias group. Two different values inside one record are
   % also a conflict.

   arguments
      existing (1, 1) struct
      incoming (1, 1) struct
   end

   % Compare each field on its own. Descriptive family and product text is not
   % the same value as a versioned source or product identifier.
   tf = false;
   fields = ["kind", "family", "source_family", "source", "source_id", ...
      "station", "product", "product_id", "relationship", "doi", ...
      "bundle_doi", "schema", "schema_version"];
   for field = fields
      [old_value, old_known] = scalarIdentity(existing, field);
      [new_value, new_known] = scalarIdentity(incoming, field);
      if old_known && new_known && old_value ~= new_value
         return
      end
   end

   % Collapse only the sampling-method alias pair that production and repaired
   % artifacts use. All other identity fields stay independent.
   [old_value, old_known, old_valid] = aliasIdentity( ...
      existing, ["sample_method", "method"]);
   [new_value, new_known, new_valid] = aliasIdentity( ...
      incoming, ["sample_method", "method"]);
   if ~old_valid || ~new_valid ...
         || (old_known && new_known && old_value ~= new_value)
      return
   end
   tf = true;
end

function [value, known, valid] = aliasIdentity(metadata, fields)
   %ALIASIDENTITY Collapse one documented field-alias group to a scalar value.
   values = strings(1, numel(fields));
   n_values = 0;
   for field = fields
      [candidate, present] = scalarIdentity(metadata, field);
      if present
         n_values = n_values + 1;
         values(n_values) = candidate;
      end
   end
   values = unique(values(1:n_values), 'stable');
   known = ~isempty(values);
   valid = numel(values) <= 1;
   value = "";
   if valid && known
      value = values(1);
   end
end

function [value, known] = scalarIdentity(metadata, field)
   %SCALARIDENTITY Normalize one concrete scalar identity value to text.
   value = "";
   known = false;
   name = char(field);
   if ~isfield(metadata, name)
      return
   end
   raw = metadata.(name);
   if (isnumeric(raw) || islogical(raw)) && isscalar(raw) && isfinite(double(raw))
      value = string(raw);
      known = true;
   elseif (ischar(raw) || isstring(raw)) && isscalar(string(raw))
      value = strtrim(string(raw));
      known = strlength(value) > 0;
   end
end
