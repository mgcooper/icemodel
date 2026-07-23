function command = fixtureFetchCommand(version, capabilities, kwargs)
   %FIXTUREFETCHCOMMAND Build a copy-paste release-data repair command.
   %
   %  command = icemodel.verification.setup.fixtureFetchCommand(version, ...
   %     capabilities, root=..., manifest=..., release_url=..., repo=...)
   %
   % Explicit nondefault roots, manifests, and release sources remain in the
   % generated command so provisioning repeats the same inspected request.

   arguments
      version (1, 1) string
      capabilities string
      kwargs.root (1, 1) string = defaultTestDataRoot()
      kwargs.manifest (1, 1) string = defaultManifestFile()
      kwargs.release_url (1, 1) string = ""
      kwargs.repo (1, 1) string = "mgcooper/icemodel"
   end

   % Keep capability order stable and quote every MATLAB string literal.
   quoted_capabilities = matlabStringLiteral( ...
      reshape(capabilities, 1, []));
   options = "capabilities=[" + strjoin(quoted_capabilities, ",") + "]";

   % Preserve only nondefault overrides; spelling out defaults adds no repair
   % information and would make the normal public command machine-specific.
   if kwargs.root ~= defaultTestDataRoot()
      options(end + 1) = "root=" + matlabStringLiteral(kwargs.root);
   end
   if kwargs.manifest ~= defaultManifestFile()
      options(end + 1) = "manifest=" ...
         + matlabStringLiteral(kwargs.manifest);
   end
   if strlength(kwargs.release_url) > 0
      options(end + 1) = "release_url=" ...
         + matlabStringLiteral(kwargs.release_url);
   end
   if kwargs.repo ~= "mgcooper/icemodel"
      options(end + 1) = "repo=" + matlabStringLiteral(kwargs.repo);
   end
   options(end + 1) = "download=true";

   command = "icemodel.verification.setup.fetchFixtures(" ...
      + matlabStringLiteral(version) + ", " + strjoin(options, ", ") + ")";
end

function literal = matlabStringLiteral(value)
   %MATLABSTRINGLITERAL Quote strings and double embedded quote characters.
   quote = string(char(34));
   literal = quote + replace(string(value), quote, quote + quote) + quote;
end

function pathname = defaultTestDataRoot()
   %DEFAULTTESTDATAROOT Canonical release-provisioned test data root.
   pathname = string(icemodel.internal.fullpath('test', 'data'));
end

function pathname = defaultManifestFile()
   %DEFAULTMANIFESTFILE Tracked authoritative release-data manifest.
   pathname = string(icemodel.internal.fullpath('test', 'assets', ...
      'icemodel-v1.1-data-manifest.json'));
end
