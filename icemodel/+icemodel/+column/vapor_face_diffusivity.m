function De_faces = vapor_face_diffusivity(De, De_N, fn)
   %VAPOR_FACE_DIFFUSIVITY Interpolate vapor diffusivity onto the faces.
   %
   %  De_faces = icemodel.column.vapor_face_diffusivity(De, De_N, fn)
   %
   % Returns the fn-weighted harmonic mean of the node diffusivities at each
   % face, Patankar (1980) Eq. 4.9. The harmonic mean makes a face between a
   % conducting cell and an insulating one behave like the insulating one. An
   % arithmetic mean does not.
   %
   % No porosity factor multiplies De. SNTHERM's De is a whole-medium
   % effective diffusivity and its heat equation carries no f_air on the
   % vapor term (Jordan 1991 Eqs. 20-21).
   %
   % One function owns this rule. The vapor mass flux and the vapor energy
   % conductance both call it. The discrete identity energy = L * mass holds
   % only while they use the same face diffusivity.
   %
   % The north and south boundary nodes follow the convention in
   % `assemble_enthalpy_system`: the caller supplies the north value, and the
   % south face repeats the deepest node, which closes it.
   %
   % Inputs
   %   De   - Diffusivity at the column nodes [m2 s-1] (JJ x 1).
   %   De_N - Diffusivity at the north boundary node [m2 s-1]. The mass flux
   %          passes the surface value; the energy side passes De(1), which
   %          zeroes that face anyway.
   %   fn   - Interface interpolation weights (JJ+1 x 1).
   %
   % Output
   %   De_faces - Diffusivity at the faces [m2 s-1] (JJ+1 x 1).
   %
   % See also: icemodel.column.vapor_face_quantities,
   %  icemodel.column.vapor_face_conductance
   %
   %#codegen

   S = numel(De);
   De_faces = 1.0 ./ ((1.0 - fn) ./ [De_N; De] + fn ./ [De; De(S)]);
end
