function link = fixtureCallerSymlink(pathname)
   %FIXTURECALLERSYMLINK Find the first caller-controlled link in a path.
   %
   %  link = icemodel.verification.setup.fixtureCallerSymlink(pathname)
   %
   % The lexical path is inspected before canonicalization so a caller-created
   % ancestor cannot redirect fixture reads or writes. Root-level macOS aliases
   % into /private are operating-system paths, not caller-controlled escapes,
   % and remain valid.

   arguments
      pathname (1, 1) string
   end

   % Normalize dot components and relative paths without resolving links.
   empty_names = javaArray('java.lang.String', 0);
   path_object = java.nio.file.Paths.get(char(pathname), empty_names);
   path_object = path_object.toAbsolutePath().normalize();
   cursor = path_object.getRoot();
   link = "";

   % Stop at the first real symbolic link so callers can issue their own
   % operation-specific error without following the redirected suffix.
   for k = 1:path_object.getNameCount()
      cursor = cursor.resolve(path_object.getName(k - 1));
      if java.nio.file.Files.isSymbolicLink(cursor) ...
            && ~isMacSystemAlias(cursor)
         link = string(cursor.toString());
         return
      end
   end
end

function tf = isMacSystemAlias(path_object)
   %ISMACSYSTEMALIAS Recognize root aliases such as /tmp -> /private/tmp.
   tf = false;
   if ~ismac || ~path_object.getParent().equals(path_object.getRoot())
      return
   end

   % Derive the trusted target from the alias name instead of maintaining a
   % platform-version-specific list of macOS aliases.
   leaf = string(path_object.getFileName().toString());
   expected = fullfile(string(path_object.getRoot().toString()), ...
      "private", leaf);
   canonical = string(java.io.File( ...
      char(path_object.toString())).getCanonicalPath());
   tf = canonical == expected;
end
