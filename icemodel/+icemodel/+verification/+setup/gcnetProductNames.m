function products = gcnetProductNames()
   %GCNETPRODUCTNAMES Return canonical Vandecrux/GC-Net product selectors.
   %
   %  products = icemodel.verification.setup.gcnetProductNames()
   %
   % Fetch validation, inventory discovery, and gcnetProductSpec share this
   % ordered registry so accepted selectors cannot drift from DOI/file metadata.

   products = ["surface", "firn_temperature", "simulated_firn"];
end
