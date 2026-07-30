function [source_dir, station] = gcnetVandecruxInputs(source_dir, station)
   %GCNETVANDECRUXINPUTS Resolve shared Vandecrux/GC-Net roots and aliases.
   source_dir = icemodel.forcing.helpers.verificationSourceDir( ...
      source_dir, "gcnet");

   station = icemodel.forcing.helpers.gcnetVandecruxStation(station);
end
