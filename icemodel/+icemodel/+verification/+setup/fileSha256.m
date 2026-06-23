function digest = fileSha256(pathname)
   %FILESHA256 Return the lowercase hex SHA-256 of a file's bytes.
   %
   %  digest = icemodel.verification.setup.fileSha256(pathname)
   %
   %  Content hash used to make the fixture bundle verifiable: packFixtures
   %  records each fixture file's SHA-256 in the bundle manifest, and
   %  fetchFixtures re-hashes the on-disk (committed or extracted) file and
   %  compares, so a corrupt or stale fixture is detected rather than silently
   %  trusted.
   %
   %  Uses java.security.MessageDigest, which ships with every MATLAB JVM, so no
   %  toolbox or shell dependency is introduced (Simulink.getFileChecksum would
   %  pull in Simulink; a shasum shell-out would be platform-fragile in CI).
   %
   %  Input
   %    pathname : string
   %        Path to an existing file.
   %
   %  Returns
   %    digest : string
   %        64-character lowercase hex SHA-256 digest.
   %
   % See also: icemodel.verification.setup.packFixtures,
   %  icemodel.verification.setup.fetchFixtures

   arguments
      pathname (1, 1) string
   end

   % Read the whole file as raw bytes; uint8 so the digest is over content, not
   % a text decoding that could vary by encoding.
   fid = fopen(pathname, 'r');
   if fid < 0
      error('icemodel:verification:fileSha256:cannotOpen', ...
         'Cannot open file for hashing: %s', pathname);
   end
   cleaner = onCleanup(@() fclose(fid));
   bytes = fread(fid, Inf, '*uint8');

   % Hash via the JVM's SHA-256 and format the signed Java byte array as
   % unsigned lowercase hex.
   md = java.security.MessageDigest.getInstance('SHA-256');
   raw = typecast(md.digest(bytes), 'uint8');
   digest = string(lower(reshape(dec2hex(raw, 2)', 1, [])));
end
