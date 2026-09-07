function assertPromiceFilledArtifact(filename, met, site, opts)
   %ASSERTPROMICEFILLEDARTIFACT Prove product and station provenance.
   %
   %  icemodel.forcing.reconstruct.assertPromiceFilledArtifact( ...
   %     filename, met, site)
   %  icemodel.forcing.reconstruct.assertPromiceFilledArtifact( ...
   %     filename, met, site, opts)
   %
   % The check is shared by runtime loading and scientific-readiness audit.
   % It rejects a stale policy digest, a gapfill_registry field that does not
   % match provenanceCodes, and incomplete per-channel provenance before a
   % promice_filled payload can be treated as model forcing.
   %
   % Note: the artifact carries gapfill_engine_version to name its producer.
   % The check requires that stamp to be a nonblank scalar and does not
   % compare it with the model version. A shipped forcing file predates the
   % release that ships it, so its engine version need not equal the current
   % model version. Requiring equality would invalidate every shipped file on
   % every version bump. Product name, policy digest, registry, and
   % per-channel provenance carry the staleness signal instead.
   %
   % Input
   %  filename  Path of the candidate promice_filled file.
   %  met       Timetable loaded from that file.
   %  site      Station token the caller is running.
   %  opts      Optional run options. This function reads only
   %            promice_filled_expected_policy_sha256, which pins a frozen
   %            release policy digest. Omit opts, or leave the field blank,
   %            to compare against the live policySha256().
   %
   % See also: icemodel.loadmet,
   %  icemodel.forcing.reconstruct.verifyPromiceFilledReadiness

   % The expected filename pattern is part of product identity; a
   % caller-supplied forcing label cannot turn an unrelated payload into
   % promice_filled.
   [~, name, extension] = fileparts(string(filename));
   site = lower(string(site));
   expected_policy_sha256 = ...
      icemodel.forcing.reconstruct.policySha256();
   if nargin >= 4 ...
         && isfield(opts, 'promice_filled_expected_policy_sha256') ...
         && ~isblanktext(opts.promice_filled_expected_policy_sha256)
      expected_policy_sha256 = string( ...
         opts.promice_filled_expected_policy_sha256);
   end
   filename_ok = startsWith(lower(string(name)), ...
      "met_" + site + "_promice_filled_") ...
      && endsWith(lower(string(name)), "_15m") ...
      && lower(string(extension)) == ".mat";

   % The filled producer stamps both identities inside the artifact; neither a
   % readiness ledger nor a caller-supplied path can substitute for them. The
   % policy digest must equal the digest registered for the selected baseline,
   % not merely a string that resembles a digest. The engine version must be
   % present and nonblank so forensics can name the producer.
   metadata = met.Properties.UserData;
   identity_fields = ["gapfill_product", "gapfill_engine_version", ...
      "gapfill_policy_sha256", "gapfill_donors", "gapfill_channels"];
   planned_channels = strings(0, 0);
   if isstruct(metadata) && isfield(metadata, 'gapfill_channels')
      planned_channels = string(metadata.gapfill_channels);
   end
   product_ok = isstruct(metadata) ...
      && all(isfield(metadata, identity_fields)) ...
      && isscalar(string(metadata.gapfill_product)) ...
      && string(metadata.gapfill_product) == "promice_filled" ...
      && isscalar(string(metadata.gapfill_engine_version)) ...
      && strlength(string(metadata.gapfill_engine_version)) > 0 ...
      && isscalar(string(metadata.gapfill_policy_sha256)) ...
      && string(metadata.gapfill_policy_sha256) ...
      == expected_policy_sha256 ...
      && isrow(planned_channels) ...
      && all(strlength(planned_channels) > 0) ...
      && numel(unique(planned_channels)) == numel(planned_channels);
   site_ok = isstruct(metadata) && isfield(metadata, 'site') ...
      && isscalar(string(metadata.site)) ...
      && lower(string(metadata.site)) == site;
   if ~site_ok
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         'filled artifact site identity does not match case %s: %s', ...
         site, filename);
   end
   if ~(filename_ok && product_ok)
      error('icemodel:loadmet:promiceFilledIdentityMismatch', ...
         ['file is not the registered promice_filled product for ' ...
         '%s: %s'], site, filename);
   end

   % Registry identity and per-sample provenance must close together; a
   % complete-looking metadata stamp cannot excuse missing channel evidence.
   if ~hasValidPromiceFilledProvenance(met, metadata)
      error('icemodel:loadmet:promiceFilledProvenanceMismatch', ...
         'file lacks complete canonical reconstruction provenance: %s', ...
         filename);
   end
end

function valid = hasValidPromiceFilledProvenance(met, metadata)
   %HASVALIDPROMICEFILLEDPROVENANCE Verify gapfill codes per channel.
   codes = icemodel.forcing.reconstruct.provenanceCodes();
   variables = string(met.Properties.VariableNames);
   channels = unique([string(metadata.gapfill_channels), ...
      icemodel.forcing.reconstruct.icemodelRequiredChannels(), ...
      icemodel.forcing.helpers.precipitationVariables()], 'stable');
   if ismember("boom_height", variables)
      channels(end + 1) = "boom_height";
   end

   % metadata.gapfill_registry must equal the append-only code list that
   % provenanceCodes returns.
   valid = isstruct(metadata) && isfield(metadata, 'gapfill_registry') ...
      && isequal(metadata.gapfill_registry, codes);
   if ~valid
      return
   end

   % Every product channel needs a uint8 code on every sample, with missing
   % reserved exactly for nonfinite values.
   allowed = struct2array(codes);
   for channel = channels
      provenance_name = channel + "_provenance";
      if ~all(ismember([channel, provenance_name], variables))
         valid = false;
         return
      end
      values = met.(channel);
      provenance = met.(provenance_name);
      missing = ~isfinite(values);
      if ~isa(provenance, 'uint8') ...
            || any(~ismember(provenance, allowed)) ...
            || any(provenance(missing) ~= codes.missing) ...
            || any(provenance(~missing) == codes.missing)
         valid = false;
         return
      end
   end
end
