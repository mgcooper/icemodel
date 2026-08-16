function [U_vap_faces, De_faces] = vapor_face_quantities( ...
      ro_vap, ro_vap_s, De, De_s, delz, fn)
   %VAPOR_FACE_QUANTITIES Diffusive vapor flux at the control-volume faces.
   %
   %  [U_vap_faces, De_faces] = icemodel.column.vapor_face_quantities( ...
   %     ro_vap, ro_vap_s, De, De_s, delz, fn)
   %
   % Solves Fick's law on the saturation vapor density at each face:
   %
   %   U_vap = -De * d(ro_vap)/dz   [kg m-2 s-1]
   %
   % Face diffusivities use the fn-weighted harmonic mean of the node values,
   % Patankar (1980) Eq. 4.9. The gradient is the secant difference of ro_vap
   % across the face, which is the conservative form: ro_vap is exponential in
   % temperature, so a node-tangent slope and a secant difference do not agree.
   %
   % No porosity factor multiplies De. SNTHERM's De is a whole-medium
   % effective diffusivity and its heat equation carries no f_air on the vapor
   % term (Jordan 1991 Eqs. 20-21).
   %
   % One function owns this face rule, so the vapor mass flux and the vapor
   % energy flux can be built from the same face quantities. That is what
   % makes the discrete identity energy = L * mass hold rather than hold
   % approximately.
   %
   % Inputs
   %   ro_vap   - Saturation vapor density at the nodes [kg m-3] (JJ x 1).
   %   ro_vap_s - Saturation vapor density at the surface [kg m-3].
   %   De       - Effective vapor diffusivity at the nodes [m2 s-1] (JJ x 1).
   %   De_s     - Effective vapor diffusivity at the surface [m2 s-1].
   %   delz     - Distance between adjacent node centers [m] (JJ+1 x 1).
   %   fn       - Interface interpolation weights (JJ+1 x 1).
   %
   % Outputs
   %   U_vap_faces - Vapor mass flux at the faces [kg m-2 s-1] (JJ+1 x 1).
   %                 Positive is downward, into the column. The bottom face
   %                 is zero: the column exchanges vapor with the surface
   %                 alone.
   %   De_faces    - Face diffusivity [m2 s-1] (JJ+1 x 1).
   %
   % References
   %   Jordan (1991), CRREL Special Report 91-16, Eqs. 20-21.
   %   Patankar (1980), Numerical Heat Transfer and Fluid Flow, Eq. 4.9.
   %
   % See also: icemodel.column.vapor_mass_transfer,
   %  icemodel.vapor.vapor_thermal_conductivity
   %
   %#codegen

   JJ = numel(ro_vap);

   % Pad with the surface above and a repeat of the deepest node below, so
   % the face arrays run over all JJ+1 interfaces.
   ro_vap_nodes = [ro_vap_s; ro_vap; ro_vap(JJ)];

   De_faces = icemodel.column.vapor_face_diffusivity(De, De_s, fn);

   U_vap_faces = -De_faces .* ...
      (ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1)) ./ delz;

   % Bottom boundary: zero flux. The padded repeat above would otherwise give
   % a zero gradient anyway; setting it makes the closed boundary explicit.
   U_vap_faces(JJ+1) = 0;
end
