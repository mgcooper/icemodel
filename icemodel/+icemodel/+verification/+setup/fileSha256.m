function digest = fileSha256(pathname)
   %FILESHA256 Return the lowercase hex SHA-256 of a file's bytes.
   %
   %  digest = icemodel.verification.setup.fileSha256(pathname)
   %
   %  Content hash that makes the fixture bundle verifiable. packFixtures
   %  records the SHA-256 of each fixture file in the bundle manifest.
   %  fetchFixtures hashes the on-disk file again, committed or extracted, and
   %  compares the two values. The check therefore detects a corrupt or stale
   %  fixture instead of trusting it.
   %
   %  Hashing itself is done by icemodel.verification.setup.bytesSha256, which
   %  uses java.security.MessageDigest. That ships with every MATLAB JVM, so
   %  this function needs no toolbox and no shell command.
   %  Simulink.getFileChecksum would require Simulink. A shasum shell command
   %  would be fragile across CI platforms.
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
   %  icemodel.verification.setup.fetchFixtures,
   %  icemodel.verification.setup.bytesSha256

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

   digest = icemodel.verification.setup.bytesSha256(bytes);
end
