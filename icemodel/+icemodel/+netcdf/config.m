function out = config()
   %CONFIG Configure icemodel.netcdf API preferences
   %
   %  CFG = ICEMODEL.NETCDF.CONFIG() Configures project preferences for the
   %  icemodel.netcdf API. Preference values are set in the 'icemodel'
   %  preference group.
   %
   % See also: ICEMODEL.NETCDF.CREATE ICEMODEL.CONFIG

   % Get environment variables
   CONTACT = icemodel.internal.contact();
   VERSION = icemodel.internal.version();
   REFERENCES = icemodel.internal.reference();

   % Create the configuration struct
   cfg.ICEMODEL_NC_TITLE = ['IceModel v' VERSION];
   cfg.ICEMODEL_NC_CONVENTIONS = 'CF-1.11';
   cfg.ICEMODEL_NC_INSTITUTION = 'UCLA';
   cfg.ICEMODEL_NC_SOURCE = ['IceModel v' VERSION];
   cfg.ICEMODEL_NC_REFERENCES = REFERENCES;
   cfg.ICEMODEL_NC_CONTACT = CONTACT;
   cfg.ICEMODEL_NC_DEFAULT_FORMAT = 'NC_FORMAT_NETCDF4';

   % Use setenv b/c setpref is way too slow
   for field = fieldnames(cfg)'
      setenv(field{:}, cfg.(field{:}))
   end

   if nargout > 0
      out = cfg;
   end
end
