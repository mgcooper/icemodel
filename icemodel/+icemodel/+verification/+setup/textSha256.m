function digest = textSha256(text)
   %TEXTSHA256 Return the lowercase hex SHA-256 of a text value's UTF-8 bytes.
   %
   %  digest = icemodel.verification.setup.textSha256(text)
   %
   %  Content hash for a value that exists in memory rather than on disk, such
   %  as the JSON encoding of a resolved options struct. fileSha256 hashes the
   %  same way, so a file digest and a text digest of identical bytes agree.
   %
   %  Hashing itself is done by icemodel.verification.setup.bytesSha256, which
   %  uses java.security.MessageDigest. That ships with every MATLAB JVM, so
   %  this function needs no toolbox and no shell command.
   %
   %  Input
   %    text : string
   %        Text to hash. Encoded as UTF-8 before hashing, so the digest does
   %        not depend on the platform's default character encoding.
   %
   %  Returns
   %    digest : string
   %        64-character lowercase hex SHA-256 digest.
   %
   % See also: icemodel.verification.setup.fileSha256,
   %  icemodel.verification.setup.bytesSha256

   arguments
      text (1, 1) string
   end

   digest = icemodel.verification.setup.bytesSha256( ...
      unicode2native(char(text), 'UTF-8'));
end
