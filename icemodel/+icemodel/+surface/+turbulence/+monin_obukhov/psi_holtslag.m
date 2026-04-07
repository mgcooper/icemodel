function psi = psi_holtslag(zeta)
   %PSI_HOLTSLAG Holtslag and de Bruin stable profile correction.
   %
   % This stable Monin-Obukhov profile correction is used for both momentum
   % and scalars in the current bulk-MO implementation.

   %#codegen

   persistent aa bb cc dd
   if isempty(aa)
      [aa, bb, cc, dd] = icemodel.parameterLookup( ...
         'thf_bulk_holtslag_aa', 'thf_bulk_holtslag_bb', ...
         'thf_bulk_holtslag_cc', 'thf_bulk_holtslag_dd');
   end

   psi = -(aa * zeta + bb * (zeta - cc / dd) * exp(-dd * zeta) ...
      + bb * cc / dd);
end
