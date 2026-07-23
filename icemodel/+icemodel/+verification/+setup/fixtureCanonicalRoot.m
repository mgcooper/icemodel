function root = fixtureCanonicalRoot(root)
   %FIXTURECANONICALROOT Return one platform-canonical absolute fixture root.
   %
   %  root = icemodel.verification.setup.fixtureCanonicalRoot(root)
   %
   % Canonicalization delegates platform aliases, separators, and traversal to
   % the filesystem. Pack/fetch separately reject symbolic links inside this
   % trusted root where ordinary file predicates would otherwise follow them.

   arguments
      root (1, 1) string
   end

   root = string(java.io.File(char(root)).getCanonicalPath());
end
