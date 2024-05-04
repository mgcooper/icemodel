% test_cvconvert

%% Conversion: Volumetric fraction to Mass
% Test values
[dz, ro_liq, ro_ice, f_liq, f_ice] = deal(1, 1000, 917, 0.5, 0.5);
[m_liq, m_ice] = icemodel.cvconvert('volumefraction', 'mass', dz, [ro_liq, ro_ice], f_liq, f_ice);
assert(isequal(m_liq, f_liq .* ro_liq .* dz), 'Test Failed: volumefraction to mass (m_liq)');
assert(isequal(m_ice, f_ice .* ro_ice .* dz), 'Test Failed: volumefraction to mass (m_ice)');

%% Conversion: Volumetric fraction to Bulk Density
[dz, ro_liq, ro_ice, f_liq, f_ice] = deal(1, 1000, 917, 0.5, 0.5);
[g_liq, g_ice] = icemodel.cvconvert('volumefraction', 'bulkdensity', dz, [ro_liq, ro_ice], f_liq, f_ice);
assert(isequal(g_liq, f_liq .* ro_liq), 'Test Failed: volumefraction to bulkdensity (g_liq)');
assert(isequal(g_ice, f_ice .* ro_ice), 'Test Failed: volumefraction to bulkdensity (g_ice)');

%% Conversion: Bulk Density to Total Density
[dz, g_ice, g_liq, g_air] = deal(1, 1, 2, 3);  % Test values
totalDensity = icemodel.cvconvert('bulkdensity', 'totaldensity', dz, 1, g_ice, g_liq, g_air);
assert(isequal(totalDensity, g_ice + g_liq + g_air), 'Test Failed: bulkdensity to totaldensity');

%% Conversion: Bulk Density to Volumetric Fraction
[dz, ro_liq, ro_ice, g_liq, g_ice] = deal(1, 1000, 917, 500, 500);
[f_liq, f_ice] = icemodel.cvconvert('bulkdensity', 'volumefraction', dz, [ro_liq, ro_ice], g_liq, g_ice);
assert(isequal(f_liq, g_liq ./ ro_liq), 'Test Failed: bulkdensity to volumefraction (f_liq)');
assert(isequal(f_ice, g_ice ./ ro_ice), 'Test Failed: bulkdensity to volumefraction (f_ice)');

%% Conversion: Bulk Density to Mass
[dz, ro_liq, ro_ice, g_liq, g_ice] = deal(1, 1000, 917, 500, 500);
[m_liq, m_ice] = icemodel.cvconvert('bulkdensity', 'mass', dz, [ro_liq, ro_ice], g_liq, g_ice);
assert(isequal(m_liq, g_liq .* dz), 'Test Failed: bulkdensity to mass (m_liq)');
assert(isequal(m_ice, g_ice .* dz), 'Test Failed: bulkdensity to mass (m_ice)');

%% Conversion: Mass to Volumetric Fraction
[dz, ro_liq, ro_ice, m_liq, m_ice] = deal(1, 1000, 917, 100, 100);
[f_liq, f_ice] = icemodel.cvconvert('mass', 'volumefraction', dz, [ro_liq, ro_ice], m_liq, m_ice);
assert(isequal(f_liq, m_liq ./ (ro_liq .* dz)), 'Test Failed: mass to volumefraction (f_liq)');
assert(isequal(f_ice, m_ice ./ (ro_ice .* dz)), 'Test Failed: mass to volumefraction (f_ice)');

%% Conversion: Mass to Bulk Density
[dz, ro_liq, ro_ice, m_liq, m_ice] = deal(1, 1000, 917, 100, 100);
[g_liq, g_ice] = icemodel.cvconvert('mass', 'bulkdensity', dz, [ro_liq, ro_ice], m_liq, m_ice);
assert(isequal(g_liq, m_liq ./ dz), 'Test Failed: mass to bulkdensity (g_liq)');
assert(isequal(g_ice, m_ice ./ dz), 'Test Failed: mass to bulkdensity (g_ice)');

%% Conversion: Volumetric fraction to Volume
[dz, f_liq, f_ice] = deal(1, 0.5, 0.5);
[v_liq, v_ice] = icemodel.cvconvert('volumefraction', 'volume', dz, [], f_liq, f_ice);
assert(isequal(v_liq, f_liq .* dz), 'Test Failed: volumefraction to volume (v_liq)');
assert(isequal(v_ice, f_ice .* dz), 'Test Failed: volumefraction to volume (v_ice)');

%% Conversion: Bulk Density to Volume
[dz, ro_liq, ro_ice, g_liq, g_ice] = deal(1, 1000, 917, 500, 500);
[v_liq, v_ice] = icemodel.cvconvert('bulkdensity', 'volume', dz, [ro_liq, ro_ice], g_liq, g_ice);
assert(isequal(v_liq, g_liq .* dz ./ ro_liq), 'Test Failed: bulkdensity to volume (v_liq)');
assert(isequal(v_ice, g_ice .* dz ./ ro_ice), 'Test Failed: bulkdensity to volume (v_ice)');

%% Conversion: Mass to Volume
[dz, ro_liq, ro_ice, m_liq, m_ice] = deal(1, 1000, 917, 100, 100);
[v_liq, v_ice] = icemodel.cvconvert('mass', 'volume', dz, [ro_liq, ro_ice], m_liq, m_ice);
assert(isequal(v_liq, m_liq ./ ro_liq), 'Test Failed: mass to volume (v_liq)');
assert(isequal(v_ice, m_ice ./ ro_ice), 'Test Failed: mass to volume (v_ice)');

%% Conversion: Mass to Mass Fraction
[dz, m_liq, m_ice] = deal(1, 100, 100);
totalMass = (m_liq + m_ice);
[fm_liq, fm_ice] = icemodel.cvconvert('mass', 'massfraction', dz, [1 1], m_liq, m_ice);
assert(isequal(fm_liq, m_liq/totalMass), 'Test Failed: massfraction to mass (m_liq)');
assert(isequal(fm_ice, m_ice/totalMass), 'Test Failed: massfraction to mass (m_ice)');

% %% Conversion: Mass Fraction to Mass
% [fm_liq, fm_ice] = deal(0.3, 0.7);  % Test values
% [m_liq, m_ice] = icemodel.cvconvert('massfraction', 'mass', dz, [ro_liq, ro_ice], fm_liq, fm_ice);
% assert(isequal(m_liq, fm_liq .* dz .* ro_liq), 'Test Failed: massfraction to mass (m_liq)');
% assert(isequal(m_ice, fm_ice .* dz .* ro_ice), 'Test Failed: massfraction to mass (m_ice)');
