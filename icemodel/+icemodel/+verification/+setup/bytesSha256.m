function digest = bytesSha256(bytes)
   %BYTESSHA256 Return the lowercase hex SHA-256 of a byte vector.
   %
   %  digest = icemodel.verification.setup.bytesSha256(bytes)
   %
   %  One hashing implementation shared by the file and text digests, so a
   %  file digest and an in-memory digest of identical bytes agree.
   %
   %  Uses java.security.MessageDigest, which ships with every MATLAB JVM, so
   %  this function needs no toolbox and no shell command.
   %
   %  Input
   %    bytes : uint8
   %        Bytes to hash, in any array shape. A file read and a text
   %        encoding return different orientations, so this function accepts
   %        both and linearizes.
   %
   %  Returns
   %    digest : string
   %        64-character lowercase hex SHA-256 digest.
   %
   % See also: icemodel.verification.setup.fileSha256,
   %  icemodel.verification.setup.textSha256

   arguments
      bytes uint8
   end

   % Hash via the JVM's SHA-256 and format the signed Java byte array as
   % unsigned lowercase hex. An empty MATLAB array marshals to a Java null
   % and throws, so the no-argument digest serves the empty input. It returns
   % the digest of nothing, which is the value an empty file must produce.
   md = java.security.MessageDigest.getInstance('SHA-256');
   if isempty(bytes)
      raw = typecast(md.digest(), 'uint8');
   else
      raw = typecast(md.digest(bytes(:)), 'uint8');
   end
   digest = string(lower(reshape(dec2hex(raw, 2)', 1, [])));
end
