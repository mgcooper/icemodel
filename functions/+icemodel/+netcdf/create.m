function varargout = create(filename, createMode, opts, ncprops, ncatts)
   %CREATE Create a new NetCDF file with specified properties and global attributes.
   %
   %  NCID = ICEMODEL.NETCDF.CREATE(FILENAME, OPTS, NCPROPS)
   %
   % Description
   %
   %  The function creates a new NetCDF file based on the specified format, fill
   %  mode, and global attributes. It handles the file backup if it already
   %  exists and applies best practices for NetCDF creation.
   %
   % Input Arguments
   %
   %  filename - Name of the NetCDF file to create.
   %  opts     - Structure with options affecting file creation behavior:
   %     opts.make_backups - Logical flag to backup the file if it exists.
   %
   %  ncprops  - Structure with properties and global attributes for the NetCDF file:
   %     createMode      - Format of the NetCDF file.
   %     ncprops.fillMode    - Fill mode for the NetCDF file, controlling prefill behavior.
   %     ncprops.comment     - User-supplied comment to include in the file's global attributes.
   %     Plus preferences set by icemodel.netcdf.config for title, Conventions, institution, source, references, and contact.
   %
   % Output Arguments
   %
   %  ncid - Numeric identifier of the newly created NetCDF file.
   %
   % Examples
   %
   %  ncid = icemodel.netcdf.create('example.nc', ...
   %     'make_backups', true, ...
   %     'format', 'NETCDF4', ...
   %     'fillMode', 'NC_NOFILL', ...
   %     'comment', 'This is a test file');
   %
   %  % close the file and examine the contents
   %  netcdf.close(ncid)
   %  info = ncinfo('example.nc');
   %
   %
   %  The CF Conventions state the following regarding global attributes: "A
   %  general description of a file's contents should be contained in the
   %  following attributes: title, history, institution, source, comment and
   %  references".
   %
   % See also: NETCDF.CREATE, ICEMODEL.NETCDF.CONFIG

   arguments
      filename (1, :) char

      % Note: this follows netcdf.create syntax, but netcdf.setDefaultFormat
      % is set in makencfile which negates the need to specify this, and could
      % lead to conflicting settings. However, it is not that simple, because I
      % think they are sometimes exclusive and sometimes not e.g. default format
      % may be NC_FORMAT_NETCDF4 but createMode could be 'NOCLOBBER', see the
      % bitwise OR options to specify two cmodes. For now, I am using NETCDF4
      % which is consistent with default NC_FORMAT_NETCDF4.
      createMode (1, :) char {mustBeMember(createMode, ...
         {'NOCLOBBER', 'CLOBBER', 'SHARE', '64BIT_OFFSET', ...
         'NETCDF4', 'CLASSIC_MODEL'})} = 'NETCDF4'

      % Local function options
      opts.make_backups (1, 1) logical = true

      % Netcdf property values
      ncprops.fillMode (1, :) char {mustBeMember(ncprops.fillMode, ...
         {'NC_FILL', 'NC_NOFILL'})} = 'NC_NOFILL'

      % Netcdf attributes (set defaults in icemodel.netcdf.config)
      ncatts.title (1, :) char = getenv('ICEMODEL_NC_TITLE')
      ncatts.Conventions (1, :) char = getenv('ICEMODEL_NC_CONVENTIONS')
      ncatts.institution (1, :) char = getenv('ICEMODEL_NC_INSTITUTION')
      ncatts.source (1, :) char = getenv('ICEMODEL_NC_SOURCE')
      ncatts.references (1, :) char = getenv('ICEMODEL_NC_REFERENCES')
      ncatts.contact (1, :) char = getenv('ICEMODEL_NC_CONTACT')
      ncatts.comment (1, :) char = ''
   end

   % Set the file creation history.
   ncatts.history = ['Created: ' ...
      char(datetime("now", "Format", "dd-MMM-uuuu hh:mm:ss", "TimeZone", "UTC"))];

   % Append any user-supplied comments to the default one.
   if ~isempty(ncatts.comment)
      ncatts.comment = [ncatts.comment newline 'created by: ' getenv('USER')];
   else
      ncatts.comment = ['created by: ' getenv('USER')];
   end

   % Back up the file if it exists.
   makebackup(filename, opts.make_backups);

   % Create the netcdf file with specified format.
   ncid = netcdf.create(filename, createMode);

   % Set the fill mode. Use NC_NOFILL to avoid duplicate writes of
   % fill values that are later replaced by data.
   netcdf.setFill(ncid, ncprops.fillMode);

   % Set global attributes.
   varid = netcdf.getConstant('GLOBAL');
   fields = string(fieldnames(ncatts));
   for thisfield = fields(:)'
      netcdf.putAtt(ncid, varid, thisfield, ncatts.(thisfield));
   end

   % Parse outputs
   switch nargout
      case 0
         % If ncid is not requested, close the file
         netcdf.close(ncid)

      case 1
         varargout{1} = ncid;

      case 2
         % Placeholder for chunksize syntax
   end
end

function makebackup(filename, make_backups)
   if isfile(filename)
      backupfile(filename, make_backups);
      fprintf('Moving to recycle bin: %s \n', filename)
      fprintf('Creating new file: %s \n', filename)
      status = recycle;
      recycle("on");
      delete(filename)
      recycle(status);
   end
end
