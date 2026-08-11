function root = fixtureCanonicalRoot(root)
   %FIXTURECANONICALROOT Return one platform-canonical absolute fixture root.
   %
   %  root = icemodel.verification.setup.fixtureCanonicalRoot(root)
   %
   % This function lets the filesystem resolve platform aliases, separators,
   % and traversal. Pack and fetch separately reject symbolic links inside this
   % trusted root, because the ordinary file tests would follow them.

   arguments
      root (1, 1) string
   end

   root = string(java.io.File(char(root)).getCanonicalPath());
end
