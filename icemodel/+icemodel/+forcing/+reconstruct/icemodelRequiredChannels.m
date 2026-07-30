function channels = icemodelRequiredChannels()
   %ICEMODELREQUIREDCHANNELS The POLICY A5 seven-channel icemodel set.
   %
   %  channels = icemodel.forcing.reconstruct.icemodelRequiredChannels()
   %
   % Role
   %  Single source of the A5 ready_icemodel channel set, consumed by the
   %  reconstruction ledger default (reconstruct.setopts
   %  required_channels) and the runtime forcing gate
   %  (icemodel.verifyPromiceFilledReadiness) so the two can never
   %  silently diverge. swu is derived and never required (A5/A16/B10);
   %  rainf is never required (rain inclusion is a model option, D-0b);
   %  snowfall input (ppt OR snowf) is graded separately by the
   %  ready_snowmodel verdict, never here.
   %
   % Returns
   %  channels : 1x7 string, the A5 icemodel forcing channels.
   %
   % See also: icemodel.forcing.reconstruct.setopts,
   %  icemodel.verifyPromiceFilledReadiness

   % Keep this explicit so an A5 change is a deliberate policy edit made
   % exactly once.
   channels = ["tair", "rh", "wspd", "psfc", "swd", "lwd", "albedo"];
end
