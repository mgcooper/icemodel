function mustBeCapHours(value)
   %MUSTBECAPHOURS Require a gap cap within any approved channel ceiling.

   caps = icemodel.forcing.reconstruct.interpolationCapHours();
   mustBeLessThanOrEqual(value, max(cell2mat(struct2cell(caps))));
end
