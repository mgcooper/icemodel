function volumeThickness = volFracToThickness(volumeFractions,layerThicknesses)
%VOLFRACTOTHICKNESS convert volumetric fractions of a layered medium to total
%volumetric thickness. if the data is 2-dimensional, the cumulative sum is
%returned which assumes the data represent changes in volume fractions over time
% 
%  volumeThickness = volFracToThickness(volumeFractions,layerThicknesses)
%  multiplies the volumetric fractions of constituent i in volumeFraction by the
%  total layer thicknesses in layerThicknesses and sums them to get the total
%  volumetric thickness of constituent i, returned as volumeThickness. 
% 
%  numel(layerThicknesses) must equal size(volumeFractions,1), meaning this
%  function assumes layerThicknesses is a 1-d vector and volumeFractions is a
%  n-d array with n<3 where the first dimension has the same number of elements
%  as layerThicknesses. Common example would be a vector of layer thicknesses
%  and an array of liquid water fractions where the rows represent layers and
%  the columns represent time steps. Multiplying the constituent layer fraction
%  of each layer by the total layer thickness yields the constituent layer
%  thickness, for example the liquid water thickness. Summing the layers yields
%  the constituent total volume thickness. Applying cumsum yields the cumulative
%  sum of constituent total volume thickness, for example the cumulative liquid
%  water thickness of a melting ice volume, or snowpack.
%  
% 
% See also 

volumeThickness = cumsum( sum( layerThicknesses(:).*volumeFractions ) );