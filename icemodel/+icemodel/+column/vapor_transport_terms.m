function [k_eff_faces, k_vap_faces, q_vap_deferred_faces, U_vap_faces, ...
      L_vap_faces] = vapor_transport_terms(T_ice, f_ice, f_liq, k_eff, ...
      ro_vap, dro_vapdT, De, delz, fn, f_res_por)
   %VAPOR_TRANSPORT_TERMS Build the coupled vapor face transport terms.
   %
   %  [k_eff_faces, k_vap_faces, q_vap_deferred_faces, U_vap_faces, ...
   %     L_vap_faces] = icemodel.column.vapor_transport_terms(T_ice, f_ice, ...
   %     f_liq, k_eff, ro_vap, dro_vapdT, De, delz, fn, f_res_por)
   %
   % Returns the vapor-coupled face terms for the enthalpy solve: the combined
   % matrix conductivity, its vapor share, the deferred vapor-energy flux, the
   % conjugate vapor-mass flux, and the face donor latent heat. The
   % conductive and diffusive node values use the same fn-weighted harmonic
   % face interpolation (Patankar 1980, Eq. 4.9).
   %
   % The exact vapor energy flux, positive downward, is
   %
   %   Q_vap = L_face * U_vap
   %         = -L_face * De_face * d(ro_vap)/dz.
   %
   % A positive donor-tangent conductivity carries the matrix part. The
   % deferred correction carries the difference between that term and the
   % exact flux. Their sum therefore stays conjugate to U_vap even across a
   % wet/dry face or an isothermal phase boundary, without putting a negative
   % coefficient in the matrix.
   %
   % K_EFF is the vapor-free node conductivity. L_face and its tangent come from
   % the donor node. The donor phase uses vapor_exchange_is_wet with the solve
   % state, and the returned L_VAP_FACES records that decision for the applied
   % mass transfer step (icemodel.column.couple_vapor_step), so the phase of the
   % transferred mass (ice or liq) matches the latent heat used for the face
   % heat flux (recomputing the phase decision after surface exchange could
   % switch the top cell's state). For the interior vapor transport computed
   % here, both boundary faces are closed; the surface exchange is computed
   % using the surface energy and mass balance instead.
   %
   % Inputs
   %   T_ice     - Node temperature [K] (JJ x 1).
   %   f_ice     - Node ice fraction [-] (JJ x 1).
   %   f_liq     - Node liquid fraction [-] (JJ x 1).
   %   k_eff     - Vapor-free conductivity at nodes [W m-1 K-1] (JJ x 1).
   %   ro_vap    - Saturation vapor density [kg m-3] (JJ x 1).
   %   dro_vapdT - Vapor-density tangent [kg m-3 K-1] (JJ x 1).
   %   De        - Effective vapor diffusivity [m2 s-1] (JJ x 1).
   %   delz      - Distance between node centers [m] (JJ+1 x 1).
   %   fn        - Interface interpolation weights [-] (JJ+1 x 1).
   %   f_res_por - Residual liquid fraction per pore volume [-].
   %
   % Outputs
   %   k_eff_faces          - Combined matrix conductivity [W m-1 K-1].
   %   k_vap_faces          - Vapor matrix conductivity [W m-1 K-1].
   %   q_vap_deferred_faces - Deferred vapor-energy flux [W m-2].
   %   U_vap_faces          - Vapor-mass flux [kg m-2 s-1].
   %   L_vap_faces          - Face donor latent heat [J kg-1] (JJ+1 x 1),
   %                          Lv where the donor node is wet, Ls otherwise.
   %
   % See also: icemodel.column.solve_column_enthalpy,
   % icemodel.column.assemble_enthalpy_system,
   % icemodel.column.vapor_exchange_is_wet
   %
   %#codegen

   persistent Ls Lv
   if isempty(Ls)
      [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   end

   JJ = numel(T_ice);

   % Pad both sides so every face follows the north/south convention used by
   % assemble_enthalpy_system. Repeating the endpoint closes its gradient.
   T_ice_nodes = [T_ice(1); T_ice; T_ice(JJ)];
   ro_vap_nodes = [ro_vap(1); ro_vap; ro_vap(JJ)];
   dro_vapdT_nodes = [dro_vapdT(1); dro_vapdT; dro_vapdT(JJ)];
   f_ice_nodes = [f_ice(1); f_ice; f_ice(JJ)];
   f_liq_nodes = [f_liq(1); f_liq; f_liq(JJ)];

   % Interpolate each node property once at the faces. The same
   % face-construction rule is used for both vapor mass and energy transport.
   k_eff_faces = 1.0 ./ ((1.0 - fn) ./ [k_eff(1); k_eff] ...
      + fn ./ [k_eff; k_eff(JJ)]);
   De_faces = 1.0 ./ ((1.0 - fn) ./ [De(1); De] ...
      + fn ./ [De; De(JJ)]);

   % Compute vapor mass flux at the cell faces [kg m-2 s-1].
   d_ro_vap = ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1);
   d_T = T_ice_nodes(2:JJ+2) - T_ice_nodes(1:JJ+1);
   U_vap_faces = -De_faces .* d_ro_vap ./ delz;

   % Select the latent heat from the node that supplies each face flux.
   wet_nodes = icemodel.column.vapor_exchange_is_wet( ...
      f_ice_nodes, f_liq_nodes, f_res_por);
   L_nodes = Ls * ones(JJ + 2, 1);
   L_nodes(wet_nodes) = Lv;
   L_north = L_nodes(1:JJ+1);
   L_south = L_nodes(2:JJ+2);
   donor_is_north = d_ro_vap <= 0;
   L_faces = L_south;
   L_faces(donor_is_north) = L_north(donor_is_north);
   L_vap_faces = L_faces;

   % Use the donor tangent for the positive matrix term.
   tangent_north = dro_vapdT_nodes(1:JJ+1);
   tangent_south = dro_vapdT_nodes(2:JJ+2);
   tangent = tangent_south;
   tangent(donor_is_north) = tangent_north(donor_is_north);
   k_vap_faces = L_faces .* De_faces .* tangent;

   % Restore the exact vapor-energy flux at the incoming Picard iterate.
   q_vap_deferred_faces = ...
      (k_vap_faces .* d_T - L_faces .* De_faces .* d_ro_vap) ./ delz;

   % Close both vapor boundaries. The ordinary conductive matrix remains active
   % at both boundaries through k_cond_faces.
   k_vap_faces(1) = 0;
   k_vap_faces(JJ+1) = 0;
   q_vap_deferred_faces(1) = 0;
   q_vap_deferred_faces(JJ+1) = 0;
   U_vap_faces(1) = 0;
   U_vap_faces(JJ+1) = 0;

   % Send the combined face conductivity to the assembler.
   k_eff_faces = k_eff_faces + k_vap_faces;
end
