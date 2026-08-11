function products = gcnetProductNames()
   %GCNETPRODUCTNAMES Return canonical Vandecrux/GC-Net product selectors.
   %
   %  products = icemodel.verification.setup.gcnetProductNames()
   %
   % Fetch validation, inventory discovery, and gcnetProductSpec share this
   % ordered registry, so the accepted selectors match the DOI and file
   % metadata.

   products = ["surface", "firn_temperature", "simulated_firn"];
end
