function k_vap_faces = vapor_face_conductance( ...
      T, f_liq, ro_vap, dro_vapdT, De, fn)
   %VAPOR_FACE_CONDUCTANCE Vapor interface conductivity at the column faces.
   %
   %  k_vap_faces = icemodel.column.vapor_face_conductance( ...
   %     T, f_liq, ro_vap, dro_vapdT, De, fn)
   %
   % Returns the interface conductivity that carries vapor latent heat,
   % built from the same face quantities the vapor mass flux uses. That makes
   % the two discretely conjugate. The energy the solve moves across a face
   % equals the mass that crosses it times the latent heat of the phase that
   % supplied it.
   %
   % The return is an interface conductivity [W m-1 K-1], Patankar's
   % interface conductivity. The conductances are the per-area aN and aS
   % [W m-2 K-1] that icemodel.column.assemble_enthalpy_system forms after
   % dividing by delz.
   %
   % The enthalpy assembly writes every face flux as an interface
   % conductivity times a temperature difference over the node spacing. The
   % vapor energy flux is
   %
   %   Q_vap = L * U_vap = -L * De_face * (ro_vap_S - ro_vap_N) / delz
   %
   % so the interface conductivity that reproduces it is
   %
   %   k_vap = L_face * De_face * (ro_vap_S - ro_vap_N) / (T_S - T_N)
   %
   % which reproduces the flux as k_vap * (T_N - T_S) / delz. The assembly
   % multiplies the north-minus-south difference, so that orientation is the
   % one to check this against, not its negative.
   %
   % The bracketed ratio is the secant slope of ro_vap over the face. It is
   % not the node-tangent derivative the default mode uses: ro_vap is
   % exponential in temperature, so the two differ across every face. Where
   % the face spans no temperature difference the secant is undefined, and
   % the tangent at the face temperature is its limit, so this uses that.
   %
   % L_face is the donor-cell choice, DesignSpec decision 5: the latent heat
   % belongs to the cell the vapor leaves, because that cell supplies the
   % phase change. A face between a wet cell and a dry one therefore carries
   % Lv when vapor flows out of the wet cell. It carries Ls when vapor flows
   % out of the dry one.
   %
   % Inputs
   %   T         - Node temperatures [K] (JJ x 1).
   %   f_liq     - Liquid fraction at the nodes [-] (JJ x 1).
   %   ro_vap    - Saturation vapor density at the nodes [kg m-3] (JJ x 1).
   %   dro_vapdT - Its temperature derivative [kg m-3 K-1] (JJ x 1), used
   %               where a face spans no temperature difference.
   %   De        - Effective vapor diffusivity at the nodes [m2 s-1] (JJ x 1).
   %   fn        - Interface interpolation weights (JJ+1 x 1).
   %
   % Output
   %   k_vap_faces - Vapor interface conductivity at the faces [W m-1 K-1]
   %                 (JJ+1 x 1). Both boundary faces are zero, matching the
   %                 mass side: the bottom is closed, and the surface face
   %                 carries the turbulent exchange rather than a diffusive
   %                 flux. The surface energy balance already accounts for
   %                 that exchange through Qe, per DesignSpec decision 2.
   %
   % See also: icemodel.column.vapor_face_quantities,
   %  icemodel.column.assemble_enthalpy_system,
   %  icemodel.vapor.latent_enthalpy_switch
   %
   %#codegen

   JJ = numel(T);

   % Pad with the surface above and a repeat of the deepest node below, the
   % same padding vapor_face_quantities uses, so the face indices agree.
   T_nodes = [T(1); T; T(JJ)];
   ro_vap_nodes = [ro_vap(1); ro_vap; ro_vap(JJ)];
   dro_vapdT_nodes = [dro_vapdT(1); dro_vapdT; dro_vapdT(JJ)];
   f_liq_nodes = [f_liq(1); f_liq; f_liq(JJ)];

   % One face diffusivity rule, shared with the mass flux. Calling the same
   % function is what makes that true. Two copies of the formula would drift
   % apart the first time one gained a weighting factor.
   De_faces = icemodel.column.vapor_face_diffusivity(De, De(1), fn);

   d_ro_vap = ro_vap_nodes(2:JJ+2) - ro_vap_nodes(1:JJ+1);
   d_T = T_nodes(2:JJ+2) - T_nodes(1:JJ+1);

   % Donor-cell latent heat. The flux runs down the vapor-density gradient,
   % so a face with more vapor below carries it upward and the lower cell is
   % the donor. Otherwise the upper cell is.
   L_nodes = icemodel.vapor.latent_enthalpy_switch(f_liq_nodes);
   L_north = L_nodes(1:JJ+1);
   L_south = L_nodes(2:JJ+2);
   donor_is_north = d_ro_vap <= 0;
   L_faces = L_south;
   L_faces(donor_is_north) = L_north(donor_is_north);

   % Secant slope where the face spans a temperature difference. Where it
   % does not, the secant is undefined and its limit is the tangent, so use
   % the donor node's derivative.
   slope = zeros(JJ+1, 1);
   spans_temperature = d_T ~= 0;
   slope(spans_temperature) = ...
      d_ro_vap(spans_temperature) ./ d_T(spans_temperature);

   dro_vapdT_north = dro_vapdT_nodes(1:JJ+1);
   dro_vapdT_south = dro_vapdT_nodes(2:JJ+2);
   tangent = dro_vapdT_south;
   tangent(donor_is_north) = dro_vapdT_north(donor_is_north);
   slope(~spans_temperature) = tangent(~spans_temperature);

   % The slope is negative at a face where the vapor flux opposes the
   % temperature gradient. A face between a colder wet cell and a warmer dry
   % one reaches that case. Ice's saturation curve sits below water's, so the
   % colder cell can hold the higher vapor density. The flux then runs up the
   % temperature gradient. That is the Bergeron-Findeisen direction, and the
   % mass flux carries it too. The conductance must therefore keep the same
   % sign to stay conjugate to the mass flux.
   %
   % The enthalpy assembly adds this to a face conduction conductivity that
   % is always positive. A negative term large enough to cancel it would cost
   % the tridiagonal system its diagonal dominance. Bead icemodel-bhk.4
   % measures how close a real column comes to that.
   k_vap_faces = L_faces .* De_faces .* slope;

   % Both boundaries carry no diffusive vapor energy. The bottom is closed.
   % The surface exchange enters through the SEB's Qe and the Neumann mass
   % flux, so a diffusive surface conductance would count it twice.
   k_vap_faces(1) = 0;
   k_vap_faces(JJ+1) = 0;
end
