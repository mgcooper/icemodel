function Pa = PRESSURE(topo)
   %PRESSURE Compute the average station pressure.
   %
   %   one_atmos  = 101300.0;
   %   scale_ht   = 8000.0;
   %   Pa         = one_atmos * exp(- topo / scale_ht);
   %
   % See also:

   Pa = 101300.0 * exp(- topo / 8000.0);
   % [Pa = J m-3 = N m-2 = kg m-1 s-2]
end
