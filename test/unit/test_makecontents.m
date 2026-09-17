function tests = test_makecontents
   %TEST_MAKECONTENTS Test the Contents.m generator of icemodel.internal.
   %
   % Each test builds a scratch namespace tree in its own temporary folder,
   % so icemodel.internal.makecontents never writes the working tree.
   tests = functiontests(localfunctions);
end

function setup(testCase)
   % A fresh scratch tree per test keeps the tests independent.
   testCase.TestData.root = string(tempname);
   mkdir(testCase.TestData.root);
   testCase.addTeardown(@() rmdir(testCase.TestData.root, 's'));
   makeFixtureTree(testCase.TestData.root);
end

function test_contents_lists_members_and_h1_lines(testCase)
   % The +pkg listing names each m-file with its dotted namespace name and
   % its H1 line, lists other files by name, and groups the rows of each
   % namespace and class folder below +pkg. A comment that holds only the
   % name, with or without a period, gives no H1 text. A comment after the
   % first code line documents that code, so it gives no H1 text either. A
   % declaration on two lines keeps the H1 line below it. A folder with no
   % file, such as +nested, gets no group. Hidden, "~", ".mex", and private
   % files get no row.

   root = testCase.TestData.root;
   icemodel.internal.makecontents(Folder=root);

   returned = readlines(fullfile(root, "+pkg", "Contents.m"), ...
      EmptyLineRule="skip");
   sep = filesep;
   expected = [
      "% +PKG"
      "%"
      "%   Contents file for +PKG and its subfolders."
      "%"
      "%   +PKG"
      "%   pkg.alpha               - First line"
      "%   pkg.beta"
      "%   notes.txt"
      "%"
      "%   +PKG" + sep + "+EMPTY"
      "%   readme.txt"
      "%"
      "%   +PKG" + sep + "+NESTED" + sep + "+DEEPER"
      "%   pkg.nested.deeper.omega - Deeper member"
      "%"
      "%   +PKG" + sep + "+SUB"
      "%   pkg.sub.delta"
      "%   pkg.sub.epsilon         - Header comment first"
      "%   pkg.sub.gamma"
      "%   pkg.sub.kappa"
      "%   pkg.sub.lambda          - Continued declaration"
      "%"
      "%   +PKG" + sep + "@THING"
      "%   pkg.thing.thing         - A class"
      "%"
      ];
   testCase.verifyEqual(returned(1:end - 1), expected);

   % The last line records when updatecontents.m wrote the file.
   testCase.verifyTrue(startsWith(returned(end), ...
      "%   updatecontents.m generated this file on "));

   % No generated line carries trailing whitespace.
   testCase.verifyFalse(any(endsWith(returned, " ")));
end

function test_contents_follow_namespace_folders_only(testCase)
   % Every namespace folder gets a Contents.m, also one with no m-file or
   % one that holds only a namespace. A private or class folder gets none.
   % A namespace folder below a plain folder takes the name of the trailing
   % "+" run.

   root = testCase.TestData.root;
   icemodel.internal.makecontents(Folder=root);

   pkg = fullfile(root, "+pkg");
   namespaces = fullfile(pkg, ["+empty", "+nested", ...
      fullfile("+nested", "+deeper"), "+sub"]);
   testCase.verifyTrue(all(arrayfun(@(f) isfile(fullfile(f, "Contents.m")), ...
      namespaces)));
   testCase.verifyFalse(isfile(fullfile(pkg, "private", "Contents.m")));
   testCase.verifyFalse(isfile(fullfile(pkg, "@thing", "Contents.m")));

   returned = readlines(fullfile(pkg, "+nested", "Contents.m"));
   testCase.verifyTrue(any(returned ...
      == "%   pkg.nested.deeper.omega - Deeper member"));

   returned = readlines(fullfile(root, "plain", "+inner", "Contents.m"));
   testCase.verifyTrue(any(returned == "%   inner.zeta - Inner member"));
end

function test_default_folder_is_the_icemodel_code_folder(testCase)
   % Without Folder, makecontents works on
   % icemodel.internal.fullpath("icemodel"). A stub fullpath on the front of
   % the path returns the scratch tree, so the test does not write the
   % working tree.

   root = testCase.TestData.root;
   stub_dir = fullfile(tempname, "stub");
   mkdir(fullfile(stub_dir, "+icemodel", "+internal"));
   testCase.addTeardown(@() rmdir(fileparts(stub_dir), 's'));
   writelines(["function p = fullpath(varargin)"
      "   p = """ + root + """;"
      "end"], fullfile(stub_dir, "+icemodel", "+internal", "fullpath.m"));
   testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stub_dir));

   icemodel.internal.makecontents();

   testCase.verifyTrue(isfile(fullfile(root, "+pkg", "Contents.m")));
end

function test_regeneration_replaces_the_listing_without_backup(testCase)
   % A later run replaces the earlier listing. The default option and
   % "-nobackup" leave no backup file behind and print nothing.

   root = testCase.TestData.root;
   icemodel.internal.makecontents(Folder=root);
   contentsfile = fullfile(root, "+pkg", "Contents.m");
   before = numel(dir(fullfile(tempdir, "pkg_Contents_*.m")));

   for call = ["icemodel.internal.makecontents(Folder=root)", ...
         "icemodel.internal.makecontents(""-nobackup"", Folder=root)"]
      writelines("% EARLIER LISTING", contentsfile, WriteMode="append");
      output = evalc(char(call));
      testCase.verifyFalse(any(readlines(contentsfile) == "% EARLIER LISTING"));
      testCase.verifyEqual(numel(dir(fullfile(tempdir, "pkg_Contents_*.m"))), ...
         before);
      testCase.verifyEqual(output, '');
   end
end

function test_backup_moves_the_earlier_listing_to_tempdir(testCase)
   % With "-backup", the earlier listing moves to a dated file in tempdir,
   % and the printed line names that file.

   root = testCase.TestData.root;
   icemodel.internal.makecontents(Folder=root);
   contentsfile = fullfile(root, "+pkg", "Contents.m");
   writelines("% EARLIER LISTING", contentsfile, WriteMode="append");
   before = tempdirBackups();

   output = evalc('icemodel.internal.makecontents("-backup", Folder=root)');

   % Every earlier listing moves to tempdir, so remove each new backup
   % before any check can stop the test. The six listings are +pkg, +empty,
   % +nested, +deeper, +sub, and +inner. The +pkg backup holds the appended
   % line.
   backups = setdiff(tempdirBackups(), before);
   testCase.addTeardown(@() arrayfun(@delete, backups));
   testCase.verifyNumElements(backups, 6);
   testCase.verifyTrue(all(arrayfun(@(b) contains(output, b), backups)));
   backup = backups(contains(backups, filesep + "pkg_Contents_"));
   testCase.assertNumElements(backup, 1);
   testCase.verifyTrue(any(readlines(backup) == "% EARLIER LISTING"));
   testCase.verifyFalse(any(readlines(contentsfile) == "% EARLIER LISTING"));
end

function test_failed_update_restores_the_earlier_listing(testCase)
   % When makecontents cannot read a member file, the update fails and the
   % error reaches the caller. A first run then leaves no Contents.m, and a later
   % run returns the earlier listing to its place.

   testCase.assumeFalse(ispc, 'The test removes read permission with chmod.');
   root = testCase.TestData.root;
   contentsfile = fullfile(root, "+pkg", "Contents.m");
   member = fullfile(root, "+pkg", "alpha.m");
   testCase.addTeardown(@() system("chmod 644 '" + member + "'"));

   % First run: no earlier listing exists, so makecontents restores nothing.
   system("chmod 000 '" + member + "'");
   testCase.verifyError(@() icemodel.internal.makecontents(Folder=root), ...
      ?MException);
   testCase.verifyFalse(isfile(contentsfile));

   % Later run: the earlier listing returns unchanged.
   system("chmod 644 '" + member + "'");
   icemodel.internal.makecontents(Folder=root);
   expected = readlines(contentsfile);
   system("chmod 000 '" + member + "'");
   testCase.verifyError(@() icemodel.internal.makecontents(Folder=root), ...
      ?MException);
   returned = readlines(contentsfile);
   testCase.verifyEqual(returned, expected);
end

function test_missing_folder_is_an_error(testCase)
   % A folder that does not exist stops the run before any file moves.

   missing = fullfile(testCase.TestData.root, "missing");
   testCase.verifyError(@() icemodel.internal.makecontents(Folder=missing), ...
      'icemodel:internal:makecontents:folderNotFound');
end

function files = tempdirBackups()
   %TEMPDIRBACKUPS Return the full paths of the Contents.m backups in tempdir.
   %
   % A backup name is <namespace>_Contents_<time>.m, so a listing before and
   % after a run finds the backups of that run.

   listing = dir(fullfile(tempdir, "*_Contents_*.m"));
   files = string(fullfile({listing.folder}, {listing.name}));
end

function makeFixtureTree(root)
   %MAKEFIXTURETREE Write a small namespace tree under ROOT.
   %
   % The tree holds one case for each listing rule: an H1 line with the
   % function name, an H1 line with no trailing period, a file with no
   % comment, a comment that is only the name, a comment that is only the
   % name and a period, a comment before the function line, a comment after
   % code, a declaration on two lines, a class folder, skipped files, a private folder, a namespace with
   % no m-file, a namespace that holds only a namespace, and a namespace
   % below a plain folder.

   pkg = fullfile(root, "+pkg");
   folders = fullfile(pkg, ["+sub", "@thing", "private", "+empty", ...
      fullfile("+nested", "+deeper")]);
   for folder = [folders, fullfile(root, "plain", "+inner")]
      mkdir(folder);
   end

   writelines(["function alpha()"; "   %ALPHA first line."; "end"], ...
      fullfile(pkg, "alpha.m"));
   writelines(["function beta()"; "end"], fullfile(pkg, "beta.m"));
   writelines("text", fullfile(pkg, "notes.txt"));
   writelines("% skipped", fullfile(pkg, "data~backup.m"));
   writelines("% skipped", fullfile(pkg, ".hidden.m"));
   writelines("binary", fullfile(pkg, "solver.mexmaca64"));
   writelines(["function delta()"; "% delta"; "end"], ...
      fullfile(pkg, "+sub", "delta.m"));
   writelines(["% header comment first."; "function epsilon()"; "end"], ...
      fullfile(pkg, "+sub", "epsilon.m"));
   writelines(["function gamma()"; "%GAMMA."; "end"], ...
      fullfile(pkg, "+sub", "gamma.m"));
   writelines(["function kappa()"; "x = 1;"; "% Internal note."; "end"], ...
      fullfile(pkg, "+sub", "kappa.m"));
   writelines(["function lambda( ..."; "   a)"; ...
      "   %LAMBDA Continued declaration."; "end"], ...
      fullfile(pkg, "+sub", "lambda.m"));
   writelines(["classdef thing"; "   %THING A class."; "end"], ...
      fullfile(pkg, "@thing", "thing.m"));
   writelines(["function helper()"; "%HELPER Private."; "end"], ...
      fullfile(pkg, "private", "helper.m"));
   writelines("text", fullfile(pkg, "+empty", "readme.txt"));
   writelines(["function omega()"; "%OMEGA Deeper member"; "end"], ...
      fullfile(pkg, "+nested", "+deeper", "omega.m"));
   writelines(["function zeta()"; "%ZETA Inner member."; "end"], ...
      fullfile(root, "plain", "+inner", "zeta.m"));
end
