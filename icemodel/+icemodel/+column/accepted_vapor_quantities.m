function [ro_vap, De] = accepted_vapor_quantities(T, f_liq)
   %ACCEPTED_VAPOR_QUANTITIES Vapor node quantities for the accepted state.
   %
   %  [ro_vap, De] = icemodel.column.accepted_vapor_quantities(T, f_liq)
   %
   % Returns the saturation vapor density and the effective vapor
   % diffusivity at the accepted substep state. The coupled vapor transport
   % uses them to build its face quantities.
   %
   % These are the two quantities icemodel.column.solve_column_enthalpy
   % evaluates on every inner iteration. The solve does not return them, so
   % this evaluates them once more at the accepted state. That is one
   % exponential and one power per accepted substep, not per iteration.
   % Returning them from the solve would remove even that, at the cost of
   % widening a signature every coupler shares. The acceptance pass decides
   % whether the cost justifies the wider signature.
   %
   % Inputs
   %   T     - Node temperatures at the accepted state [K] (JJ x 1).
   %   f_liq - Liquid fraction at the accepted state [-] (JJ x 1).
   %
   % Outputs
   %   ro_vap - Saturation vapor density [kg m-3] (JJ x 1).
   %   De     - Effective vapor diffusivity [m2 s-1] (JJ x 1).
   %
   % See also: icemodel.column.couple_vapor_transport,
   %  icemodel.vapor.saturation_vapor_density,
   %  icemodel.vapor.vapor_diffusivity
   %
   %#codegen

   ro_vap = icemodel.vapor.saturation_vapor_density(T, f_liq);
   De = icemodel.vapor.vapor_diffusivity(T);
end
