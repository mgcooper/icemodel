function [k_vap_faces, q_vap_deferred] = vapor_face_conductance( ...
      T, f_ice, f_liq, ro_vap, dro_vapdT, De, delz, fn, f_res_por)
   %VAPOR_FACE_CONDUCTANCE Vapor interface conductivity and deferred flux.
   %
   %  [k_vap_faces, q_vap_deferred] = icemodel.column.vapor_face_conductance( ...
   %     T, f_ice, f_liq, ro_vap, dro_vapdT, De, delz, fn, f_res_por)
   %
   % Returns the two pieces that carry vapor latent heat in the coupled
   % enthalpy solve, built from the same face quantities the vapor mass flux
   % uses. Their sum is what makes the energy and the mass discretely
   % conjugate: at convergence, the energy each face moves equals the mass
   % that crosses it times the latent heat of the phase that supplied it.
   %
   % The exact vapor energy flux at a face, positive downward, is
   %
   %   Q_vap = L_face * U_vap = -L_face * De_face * (ro_vap_S - ro_vap_N) / delz
   %
   % No single conductance reproduces that flux in every state. The secant
   % slope of ro_vap goes negative where the flux opposes the temperature
   % gradient: a colder wet cell can hold more vapor than a warmer dry one,
   % because ice's saturation curve sits below water's. That is the
   % Bergeron-Findeisen direction, and the mass flux carries it, so the
   % energy must carry it too. A face can also move vapor with no
   % temperature difference at all, at an isothermal wet/dry boundary,
   % where no conductance times a zero difference reproduces a nonzero
   % flux. Patankar (1980) section 7.2 splits the job:
   %
   %   1. K_VAP_FACES [W m-1 K-1] is the matrix part: the donor cell's
   %      tangent slope, L_face * De_face * dro_vapdT_donor. The tangent is
   %      positive everywhere, because saturation vapor density rises with
   %      temperature, so the assembled system keeps positive off-diagonals
   %      and its diagonal dominance at every face.
   %   2. Q_VAP_DEFERRED [W m-2], positive downward, is the deferred
   %      correction: the exact face flux minus what the matrix part moves
   %      at the incoming iterate. The assembly adds its in-minus-out
   %      difference to the source vector. The Picard loop re-evaluates it
   %      each iteration, so the converged flux is exactly L_face * U_vap
   %      while the matrix never carries the sign.
   %
   % L_face is the donor-cell choice, DesignSpec decision 5: the latent
   % heat belongs to the cell the vapor leaves, because that cell supplies
   % the phase change. The donor's phase comes from
   % icemodel.column.vapor_exchange_is_wet, the same predicate the mass
   % applier uses, so the energy and the mass never disagree about which
   % latent heat a face carries. icemodel.vapor.latent_enthalpy_switch
   % stays the owner of the storage and conduction switch inside the
   % enthalpy solve; it does not decide donor faces.
   %
   % Inputs
   %   T         - Node temperatures [K] (JJ x 1).
   %   f_ice     - Ice fraction at the nodes [-] (JJ x 1).
   %   f_liq     - Liquid fraction at the nodes [-] (JJ x 1).
   %   ro_vap    - Saturation vapor density at the nodes [kg m-3] (JJ x 1).
   %   dro_vapdT - Its temperature derivative [kg m-3 K-1] (JJ x 1).
   %   De        - Effective vapor diffusivity at the nodes [m2 s-1] (JJ x 1).
   %   delz      - Distance between adjacent node centers [m] (JJ+1 x 1).
   %   fn        - Interface interpolation weights (JJ+1 x 1).
   %   f_res_por - Residual liquid-water fraction per pore volume [-].
   %
   % Outputs
   %   k_vap_faces    - Matrix-part vapor interface conductivity at the
   %                    faces [W m-1 K-1] (JJ+1 x 1), always nonnegative.
   %   q_vap_deferred - Deferred vapor energy flux at the faces [W m-2]
   %                    (JJ+1 x 1), positive downward.
   %
   % Both boundary faces are zero in both outputs, matching the mass side:
   % the bottom is closed, and the surface face carries the turbulent
   % exchange rather than a diffusive flux. The surface energy balance
   % already accounts for that exchange through Qe, per DesignSpec
   % decision 2.
   %
   % See also: icemodel.column.vapor_face_quantities,
   %  icemodel.column.assemble_enthalpy_system,
   %  icemodel.column.vapor_exchange_is_wet
   %
   %#codegen

   persistent Ls Lv
   if isempty(Ls)
      [Ls, Lv] = icemodel.physicalConstant('Ls', 'Lv');
   end

   JJ = numel(T);

   % Pad with the surface above and a repeat of the deepest node below, the
   % same padding vapor_face_quantities uses, so the face indices agree.
   T_nodes = [T(1); T; T(JJ)];
   ro_vap_nodes = [ro_vap(1); ro_vap; ro_vap(JJ)];
   dro_vapdT_nodes = [dro_vapdT(1); dro_vapdT; dro_vapdT(JJ)];
   f_ice_nodes = [f_ice(1); f_ice; f_ice(JJ)];
   f_liq_nodes = [f_liq(1); f_liq; f_liq(JJ)];

   % One face diffusivity rule, shared with the mass flux. Calling the same
   % function is what makes that true. Two copies of the formula would drift
   % apart the first time one gained a weighting factor.
   De_faces = icemodel.column.vapor_face_diffusivity(De, De(1), fn);

   d_ro_vap = ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1);
   d_T = T_nodes(2:JJ+2) - T_nodes(1:JJ+1);

   % Donor-cell latent heat. The flux runs down the vapor-density gradient,
   % so a face with more vapor below carries it upward and the lower cell is
   % the donor. Otherwise the upper cell is. The wet/dry phase decision is
   % icemodel.column.vapor_exchange_is_wet, the mass applier's predicate,
   % so the two paths use one latent heat in the band where the residual
   % floor and the fixed storage threshold disagree.
   wet_nodes = icemodel.column.vapor_exchange_is_wet( ...
      f_ice_nodes, f_liq_nodes, f_res_por);
   L_nodes = Ls * ones(JJ + 2, 1);
   L_nodes(wet_nodes) = Lv;
   L_north = L_nodes(1:JJ+1);
   L_south = L_nodes(2:JJ+2);
   donor_is_north = d_ro_vap <= 0;
   L_faces = L_south;
   L_faces(donor_is_north) = L_north(donor_is_north);

   % Donor tangent for the matrix part. Any positive coefficient converges
   % to the same answer once the deferred term corrects it; the donor
   % tangent is the local slope of the curve the flux actually follows, so
   % it is the Newton-like choice.
   dro_vapdT_north = dro_vapdT_nodes(1:JJ+1);
   dro_vapdT_south = dro_vapdT_nodes(2:JJ+2);
   tangent = dro_vapdT_south;
   tangent(donor_is_north) = dro_vapdT_north(donor_is_north);

   % Matrix part: always nonnegative, so the assembled off-diagonals stay
   % nonnegative and the diagonal dominance the tridiagonal solve needs
   % holds at every face, in every state.
   k_vap_faces = L_faces .* De_faces .* tangent;

   % Deferred part: exact flux minus what the matrix part moves at this
   % iterate, positive downward. The matrix moves -k_vap * d_T / delz, and
   % the exact flux is -L * De * d_ro_vap / delz, so the difference below
   % restores the exact flux when the assembly adds it to the source.
   q_vap_deferred = ...
      (k_vap_faces .* d_T - L_faces .* De_faces .* d_ro_vap) ./ delz;

   % Both boundaries carry no diffusive vapor energy. The bottom is closed.
   % The surface exchange enters through the SEB's Qe and the Neumann mass
   % flux, so a diffusive surface term would count it twice.
   k_vap_faces(1) = 0;
   k_vap_faces(JJ+1) = 0;
   q_vap_deferred(1) = 0;
   q_vap_deferred(JJ+1) = 0;
end
