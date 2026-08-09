function [solid_mwe, liquid_mwe, enthalpy_j_m2] = ...
      integrate_column_budget(T, f_ice, f_liq, dz)
   %INTEGRATE_COLUMN_BUDGET Return column-integrated mass and enthalpy storage.
   %
   %  [solid_mwe, liquid_mwe] = ...
   %     icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz)
   %  [solid_mwe, liquid_mwe, enthalpy_j_m2] = ...
   %     icemodel.column.integrate_column_budget(T, f_ice, f_liq, dz)
   %
   % solid_mwe and liquid_mwe are positive storage depths in metres water
   % equivalent. enthalpy_j_m2 is the column integral of the solver's
   % bulk_enthalpy measure, using its documented dry-mixture reference and
   % omitting vapor, so remeshing is compared on one fixed material basis.
   %
   % All MWE outputs use the solver's physical intrinsic phase densities and
   % physical liquid-water density as the fixed reference. The use_ro_glc
   % option changes initialization fractions only; it does not redefine this
   % basis. See icemodel-hd6.
   %
   % DZ may be a scalar uniform-cell thickness or a vector matching the state.
   %
   %#codegen

   % Cache the constants. In the diagnostic profile this is called once per
   % accepted substep from each accumulator, plus once per merge event, and
   % physicalConstant rebuilds its whole table on every call.
   persistent ro_ice ro_liq
   if isempty(ro_ice)
      [ro_ice, ro_liq] = icemodel.physicalConstant('ro_ice', 'ro_liq');
   end

   % Integrate solid and liquid phase storage in mwe, on the same physical
   % basis as solver phase change.
   solid_mwe = ro_ice / ro_liq * sum(f_ice .* dz);
   liquid_mwe = sum(f_liq .* dz);

   % Compute the documented-reference enthalpy only when the caller needs it.
   if nargout > 2
      f_wat = icemodel.column.water_fraction(f_ice, f_liq);
      H = icemodel.column.bulk_enthalpy(T, f_ice, f_liq, f_wat);
      enthalpy_j_m2 = sum(H .* dz);
   end
end
