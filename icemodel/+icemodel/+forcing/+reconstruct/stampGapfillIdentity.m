function met = stampGapfillIdentity(met, site, codes, plan, family, donor_sites)
   %STAMPGAPFILLIDENTITY Stamp the gapfill_* identity fields on a filled met.
   %
   %  met = icemodel.forcing.reconstruct.stampGapfillIdentity( ...
   %     met, site, codes, plan, family, donor_sites)
   %
   % fillPromiceStation calls this once per station after every channel is
   % filled. assertPromiceFilledArtifact reads gapfill_product,
   % gapfill_engine_version, gapfill_policy_sha256, gapfill_donors, and
   % gapfill_channels as the artifact identity. gapfill_generated_utc records
   % when the artifact was produced; it is metadata only, so a staged
   % artifact without it still loads.
   %
   % Input
   %  met          Filled timetable whose UserData receives the stamps.
   %  site         Compact station token, such as "kanm".
   %  codes        Provenance registry from provenanceCodes.
   %  plan         Fill plan; reads plan.split.seed and plan.channels.channel.
   %  family       Product family, such as "promice".
   %  donor_sites  Donor station tokens used by the fill.
   %
   % See also: icemodel.forcing.reconstruct.fillPromiceStation,
   %  icemodel.forcing.reconstruct.assertPromiceFilledArtifact

   ud = met.Properties.UserData;
   % Runtime identity uses the compact token, not the separator-bearing
   % display name inherited from native metadata.
   ud.site = site;
   ud.gapfill_registry = codes;
   ud.gapfill_seed = plan.split.seed;
   ud.gapfill_product = char(family + "_filled");
   ud.gapfill_channels = string({plan.channels.channel});
   ud.gapfill_engine_version = string(icemodel.internal.version());
   ud.gapfill_policy_sha256 = ...
      icemodel.forcing.reconstruct.policySha256();
   ud.gapfill_donors = donor_sites(:).';
   % ISO 8601 UTC text, so the value reads the same in every locale and
   % survives a MAT round trip without a datetime time zone.
   ud.gapfill_generated_utc = string(datetime('now', 'TimeZone', 'UTC', ...
      'Format', "uuuu-MM-dd'T'HH:mm:ss'Z'"));
   met.Properties.UserData = ud;
end
