%[text] # Vapor Density
%[text] NOTE: In the expressions below I've defined $c = c - T\_f$ to simplify notation, so what is normally $c + T - T\_f$ is now $c + T$.
%%
clean

method = "buck";

syms n T T_f D_e0 e_s0 b c L_v R_v k_v n_m1 n_m2 a
syms D_e e_s rho_v d_rho_v_dT d_e_s_dT

assume([n T T_f D_e0 e_s0 b c L_v R_v k_v], 'clear')
assume([D_e e_s rho_v], 'clear')

assume([D_e0 e_s0 b c L_v R_v T_f T n], {'positive', 'real'})
assume([D_e e_s rho_v], {'positive', 'real'})
%%
%[text] Write the expressions
T_d_expr = T - T_f;
n_m1_expr = n - 1;
n_m2_expr = n - 2;
D_e_expr = D_e0 * (T / T_f) ^ n %[output:3437fefe]
if method == "buck" %[output:group:360166a3]
   % Note - it doesn't matter if symbolically c is replaced with c - Tf, b/c c
   % is only involved in this expression and is isolated from all other terms,
   % but use the c + T version for consistency with the latex notes. 
   % e_s_expr = e_s0 * exp(b * T_d_expr / (c + T_d_expr))
   e_s_expr = e_s0 * exp(b * T_d_expr / (c + T)) %[output:4f3b72fd]
elseif method == "ambaum"
   e_s_expr = a * exp(b/T) * T ^ c
end %[output:group:360166a3]
rho_v_expr = e_s_expr / (R_v * T) %[output:3c5d355f]
%%
%[text] Compute the symbolic derivatives
d_rho_v_dT_expr = diff(rho_v_expr, T, 1);
d2_rho_v_dT2_expr = diff(rho_v_expr, T, 2);
d_e_s_dT_expr = diff(e_s_expr, T, 1);
d2_e_s_dT2_expr = diff(e_s_expr, T, 2);
%%
%[text] ## Definitions
%[text] Note the following equivalent terms. Any term involving $c^\* + T$ is for buck, with $c^\* = c - T\_f$
%[text] $\\frac{d\\rho\_v}{dT} = \\frac{1}{R\_v T} \\left(\\frac{de\_s}{dT} - \\frac{e\_s}{T} \\right)$
%[text] $\\quad = \\rho\_v \\left(\\frac{1}{e\_s} \\frac{de\_s}{dT} - \\frac{1}{T} \\right)$
%[text] NOTE: I don't think the next three are correct for Buck, for sure the second and third do not follow from the first, there is an algebra error
%[text] $\\frac{d^2\\rho\_v}{dT^2} = \\frac{2 \\rho\_v}{T} \\left( \\frac{1}{T} - \\frac{1}{e\_s}\\frac{de\_s}{dT}\\right) + \\frac{\\rho\_v}{e\_s} \\frac{de\_s}{dT} \\left( \\frac{1}{e\_s} \\frac{de\_s}{dT} - \\frac{2}{c^\* + T} \\right)$
%[text] $\\quad \\quad = \\frac{2 \\rho\_v}{T} \\left( \\frac{2}{T} - \\frac{ 1 }{ \\rho\_v } \\frac{ d \\rho\_v }{ dT } \\right) + \\left( \\frac{d\\rho\_v}{dT} + \\frac{\\rho\_v}{T} \\right) \\left( \\frac{ 1 }{ \\rho\_v } \\frac{ d \\rho\_v }{ dT } + \\frac{ 1 }{ T } - \\frac{2}{c^\* + T} \\right)$
%[text] $\\quad \\quad = \\frac{2 \\rho\_v}{T} \\left( \\frac{2}{T} - \\frac{ 1 }{ \\rho\_v } \\frac{ d \\rho\_v }{ dT } \\right) - \\left( \\frac{d\\rho\_v}{dT} + \\frac{\\rho\_v}{T} \\right) \\left( \\frac{2}{c^\* + T} - \\frac{ 1 }{ T } - \\frac{ 1 }{ \\rho\_v } \\frac{ d \\rho\_v }{ dT }\\right)$
%[text] The next ones are confirmed correct and I think the first has the minimum number of calculations:
%[text] $\\quad \\quad =\\left( \\frac{\\rho\_v}{T} + \\frac{d\\rho\_v}{dT} \\right) \\left( \\frac{1}{\\rho\_v}\\frac{d\\rho\_v}{dT} - \\frac{2}{c^\*+T} \\right) - \\frac{1}{T}\\frac{d\\rho\_v}{dT} + \\frac{1}{T} \\frac{\\rho\_v}{T}$
%[text] $\\quad \\quad =\\left( \\frac{\\rho\_v}{T} + \\frac{d\\rho\_v}{dT} \\right) \\left( \\frac{1}{\\rho\_v}\\frac{d\\rho\_v}{dT} - \\frac{2}{c^\*+T} -\\frac{1}{T} \\right) + \\frac{2 \\rho\_v}{T^2}$
%[text] $\\quad \\quad = \\frac{1}{T} \\left( \\frac{\\rho\_v}{T} - \\frac{d\\rho\_v}{dT} \\right) + \\left( \\frac{\\rho\_v}{T} + \\frac{d\\rho\_v}{dT} \\right) \\left( \\frac{1}{\\rho\_v} \\frac{ d\\rho\_v }{ dT } - \\frac{ 2 }{ c^\* + T} \\right)$
%[text] $\\quad \\quad = \\frac{\\rho\_v}{T} \\left( \\frac{1}{T} - \\frac{ 2 }{ c^\* + T} \\right) + \\frac{d\\rho\_v}{dT} \\left( \\frac{1}{\\rho\_v} \\frac{ d\\rho\_v }{ dT } - \\frac{ 2 }{ c^\* + T} \\right)$
% rho_v ./ T * (1 ./ T - 2 ./ (c + T - T_f)) + d_rho_v_dT .* ( d_rho_v_dT ./ rho_v - 2 ./ (c + T - T_f))
%[text] I think above is most concise in terms of $\\rho\_v$. What is most concise in terms of $e\_s$?
%[text] This first one has the fewest number of operations and no powers. This is implemented in icemodel.vapor.vappress:
%[text] $\\frac{d^2\\rho\_v}{dT^2} = \\frac{1}{R\_v T} \\left\[ \\frac{de\_s}{dT} \\left( \\frac{1}{e\_s} \\frac{de\_s}{dT} - \\frac{2}{c^\*} \\right) - \\frac{2}{T} \\left( \\frac{de\_s}{dT} - \\frac{e\_s}{T} \\right) \\right\]$
% (de_s_dT * (de_s_dT ./ e_s - 2 ./ cstar) - 2 * (de_s_dT - es ./ T) ./ T) ./ (Rv * T)
%[text] $= \\frac{1}{R\_v T} \\left\[ \\frac{de\_s}{dT} \\left( \\frac{1}{e\_s} \\frac{de\_s}{dT} - \\frac{2}{c^\*} - \\frac{2}{T} \\right) + \\frac{2 e\_s}{T^2} \\right\]$
% (de_s_dT * (de_s_dT ./ e_s - 2 ./ cstar - 2 ./ T) + 2 * es ./ T.^2) ./ (Rv * T)
%[text] $= \\frac{1}{R\_v T} \\left\[ \\frac{de\_s}{dT} \\left( \\frac{1}{e\_s} \\frac{de\_s}{dT} - 2\\left(\\frac{1}{c^\*} + \\frac{1}{T} \\right) \\right) + \\frac{2 e\_s}{T^2} \\right\]$
% (de_s_dT * (de_s_dT ./ e_s - 2 * (1 ./ cstar + 1 ./ T)) + 2 * es ./ T.^2) ./ (Rv * T)
%[text] $= \\frac{1}{R\_v T^2} \\left\[ \\frac{de\_s}{dT} \\left( \\frac{T}{e\_s} \\frac{de\_s}{dT} - 2 \\left( \\frac{T}{c^\*} +1\\right) \\right) + \\frac{2 e\_s}{T} \\right\]$
% (de_s_dT * (de_s_dT .* T ./ e_s - 2 * (T ./ cstar + 1)) + 2 * es ./ T) ./ (Rv * T.^2)
%%
%[text] ## Useful substitution definitions
%[text] These definitions are for Buck with $c^\* = c - T\_f$.  $\\sigma$ numbering may differ from output depending on order of operations. 
%[text] $\\sigma\_1 = \\left( \\frac{ b }{ c^\* + T} - \\frac{ b\[T - T\_f\] }{ \[c^\* + T\]^2 } \\right)$
%[text] $\\quad = \\left( \\frac{ bc }{ \[c^\* + T \]^2} \\right)$
%[text] $\\quad = \\left( \\frac{ 1 }{ e\_s } \\frac{ de\_s }{ dT } \\right)$
%[text] $\\quad = \\left( \\frac{ 1 }{ \\rho\_v } \\frac{ d \\rho\_v }{ dT } + \\frac{ 1 }{ T } \\right )$
%[text] $\\sigma\_2 = \\frac{ 2b }{ \[c^\* + T\]^2 } - \\frac{ 2b\[T - T\_f\] }{ \[c^\* + T\]^3 }$
%[text] $\\quad =  \\frac{ 2bc }{ \[c^\* + T\]^3}$
%[text] $\\quad = \\left( \\frac{ 1 }{ e\_s } \\frac{ de\_s }{ dT } \\right) \\frac{ 2 }{ c^\* + T }$
%[text] $\\quad = \\left( \\frac{ 1 }{ \\rho\_v } \\frac{ d \\rho\_v }{ dT } + \\frac{ 1 }{ T } \\right)  \\frac{ 2 }{ c^\* + T }$
%[text] $\\rho\_v \\sigma\_1 = \\frac{\\rho\_v}{e\_s}\\frac{de\_s}{dT} = \\frac{d\\rho\_v}{dT} + \\frac{\\rho\_v}{T}$
%[text] These are useful for $k\_v$:
%[text] $k\_v = \\rho\_v D\_e L\_v \\left(\\sigma\_1 - \\frac{1}{T} \\right)$
%[text] $\\sigma\_1 = \\frac{k\_v}{\\rho\_v D\_e L\_v}  - \\frac{1}{T}$
%[text] $\\frac{\\rho\_v D\_e L\_v}{T} = \\frac{\\rho\_v k\_v}{T \\frac{d \\rho\_v}{dT}}$
%[text] $\\frac{k\_v}{\\rho\_v D\_e L\_v} = \\frac{1}{T} \\frac{d \\rho\_v}{dT}$
%%
%[text] ### Compute $k\_v$ and it's temperature derivative
%[text] Create $k\_v$ and its derivative dynamically from the expressions
k_v_expr = L_v * D_e_expr * d_rho_v_dT_expr;
d_k_v_dT_expr = diff(k_v_expr, T, 1);
%[text] Substitute the expressions for $D\_e$ and $\\rho\_v$. I tried subbing them individually and there is no benefit.
% It won't factor out n-1 or n-2
% d_k_v_dT = subs(d_k_v_dT_expr, {n_m1_expr, D_e_expr, rho_v_expr}, {n_m1, D_e, rho_v})
% d_k_v_dT = subs(d_k_v_dT_expr, {n_m2_expr, D_e_expr, rho_v_expr}, {n_m2, D_e, rho_v})
d_k_v_dT = subs(d_k_v_dT_expr, {D_e_expr, rho_v_expr}, {D_e, rho_v}) %[output:3d0af837]
% Sub in the alternative form for De in the first term above
d_k_v_dT = subs(d_k_v_dT, D_e0 * (T/T_f)^(n-1) / T_f, D_e / T) %[output:8722b300]
% If using Ambaum:
if method == "ambaum"
   d_k_v_dT = subs(d_k_v_dT, T^(c-1) * a * exp(b/T) / (R_v * T), rho_v/T);
   d_k_v_dT = subs(d_k_v_dT, T^(c-2) * a * exp(b/T) / (R_v * T), rho_v/T^2)
end

% Buck:
% dkdT = k_vap * (n ./ T - 2 * (ro_vap ./ (T .* drov_dT) + 1) ./ (c + T - Tf) + drov_dT ./ ro_vap + ro_vap ./ (T .* drov_dT)); % 15
% dkdT = k_vap * (n ./ T + drov_dT ./ ro_vap - 2 ./ (c + T - Tf)) + Lv * De .* ro_vap .* (1 ./ T - 2 ./ (c + T - Tf)) ./ T; % 18
% dkdT = k_vap * (n ./ T + drov_dT ./ ro_vap - 2 ./ (c + T - Tf)) + k_vap .* ro_vap .* (1 ./ T - 2 ./ (c + T - Tf)) ./ (T .* drov_dT); % 18
%%
%[text] This section experimented with subbing in (n-1) but it doesn't work as intended. 
% d_k_v_dT = subs(d_k_v_dT_expr, n_m1_expr, n_m1)
% Can't sub in D_e after nm1, but rho_v works but there's no benefit:
% d_k_v_dT = subs(d_k_v_dT, rho_v_expr, rho_v)
%%
%[text] Substitute in terms of $e\_s$
% k_v = subs(k_v_expr, e_s_expr, e_s)
% k_v = subs(k_v, D_e_expr, D_e)
%%
% Here I wanted to try getting the individual terms and then subbing them but
% the problem is the individual terms do not necessarily correspond to
% closed-form substitutions
% subtc = children(d2_rho_v_dT2_expr)
% subtc{1}
% subtc{2}
% test = subs(d2_rho_v_dT2, subtc{1}, e_s)
%%
%[text] ### Substitute variables for the symbolic expressions to simplify
%[text] Substitute in terms of $e\_s$:
d_e_s_dT = subs(d_e_s_dT_expr, e_s_expr, e_s) %[output:0fecc6e0]
%[text] Substitute in terms of $e\_s$:
% d_rho_v_dT = subs(d_rho_v_dT_expr, e_s_expr, e_s)
% d2_rho_v_dT2 = subs(d2_rho_v_dT2_expr, e_s_expr, e_s)
%[text] Substitute in terms of $\\rho\_v$:
d_rho_v_dT = subs(d_rho_v_dT_expr, rho_v_expr, rho_v) %[output:28891588]
if method == "ambaum"
   d_rho_v_dT = subs(d_rho_v_dT, T^(c-1) * a * exp(b/T) / (R_v * T), rho_v/T)
end
%%
%[text] Second derivative
% Substitute the rho_v expression
d2_rho_v_dT2 = subs(d2_rho_v_dT2_expr, rho_v_expr, rho_v) %[output:2c83f7d8]
% Substitute the alternative forms that matlab doesn't sub automatically
if method == "ambaum"
   d2_rho_v_dT2 = subs(d2_rho_v_dT2, T^(c-1) * a * exp(b/T) / (R_v * T), rho_v/T);
   d2_rho_v_dT2 = subs(d2_rho_v_dT2, T^(c-2) * a * exp(b/T) / (R_v * T), rho_v/T^2)
end
%%
%[text] ## Ambaum
%[text] Check the derivative of the Ambaum expression to see how much more complex it becomes
syms a b c e_s
e_s_expr = a * exp(b/T) * T ^ c  %[output:0de36b50]
de_s_dT_expr = diff(e_s_expr, T, 1) %[output:7e05bf15]
d2e_s_dT2_expr = diff(e_s_expr, T, 2);
%[text] Rewrite
de_s_dT = subs(de_s_dT_expr, e_s_expr, e_s);
de_s_dT = subs(de_s_dT, T^(c-1) * a * exp(b/T), e_s/T) %[output:1bd7c419]
%[text] Rewrite the second derivative in terms of the first
% Doesn't work
% syms de_s
% d2e_s_dT2 = subs(d2e_s_dT2_expr, de_s_dT_expr, de_s)
%[text] Rewrite the second derivative in terms of $e\_s$
d2e_s_dT2 = subs(d2e_s_dT2_expr, e_s_expr, e_s);
d2e_s_dT2 = subs(d2e_s_dT2, T^(c-1) * a * exp(b/T), e_s/T);
d2e_s_dT2 = subs(d2e_s_dT2, T^(c-2) * a * exp(b/T), e_s/T^2) %[output:85d99bbb]
%%
diff(exp(b/T), T) %[output:7877f4b3]
%%
% NOTE: The only functional difference between Romps and Ambaum is Romps uses
% cv_liq whereas Ambaum uses cp_liq in the exponent on (T/To). I think Ambaum
% addresses this on page 4253, top right. 

% Based on Ambaum's Fig 4, it appears fine to compute wrt water or solid, it
% will only make a large difference at low temperatures, 
% Since RH is likely relative to water, and mar definitely was computed wrt
% water, for the surface 

% Ambaum key points: 
% - Recommended values are for the entire range -60-100oC
% - Recommended values set cp_liq = 4220 and optimize cp_vap = 2040
% - Figure 3 shows less error with the actual triple point cp_vap = 1888
% - But this is supercooled water, and the difference is minor above -10oC
% - Since cp_vap is not used elsewhere, and if it were, the optimal value
%   could be adopted, there's no reason not to use cp_vap = 2040
%
% Over ice, he says "triple point values can be used":
% - cp_ice - cp_vap = 2097 - 1888.2 = 208.8
% But his recommended value for cp_ice - cp_vap is 212, so the triple point
% value of cp_vap over ice must be:
% - cp_vap = cp_ice - 212 = 2097 - 212 = 1885
% 
% 

% Ambaum Eq. 13-15 and 17-19:
% cp_liq = 4220;     % Triple-point value, this is fixed
% cp_lmv = 2180;     % cp_liq - cp_vap, optimized
% cp_vap = 2040;     % optimized; actual triple point value = 1888.2;
% cp_ice = 2097;     % Triple-point value, this is fixed
% cp_imv = 212;      % cp_ice - cp_vap = 2097 - 1885

% L0 = 2.501e6;      % Latent heat of vaporization at 0°C [J kg-1]
% es0 = 611.655;
% Rv = 461.52;       % Gas constant for water vapor [J kg-1 K-1]

% cp_

% Note: above gives cp_vap = 1879.
% Ambaum got 1888.2 at the triple point, and says this value yields a value
% for cp_liq "fairly close to its measured triple point value". NOTE: it is
% not correct to compute cp_liq = 2180 + cp_vap = 4059, because the value
% 2180 was selected by Ambaum to give cp_vap = 2040. 
% Thus cp_liq = 2180 + 2040 = 4220, almost identical to my value 4218.

% Ambaum says at low temperatures, the best fit is for the actual triple
% point value of cp_vap, not the one above, which produces better results
% over a larger temperature range up to 100oC. 
%
% Ambaum's triple-point values:
% cp_vap = 1888.2;
% cp_liq = 4220;
% cp_ice = 2097;
%
% Those values give cp_liq - cp_vap = 2331.8;
% See Ambaum Figure 3 for why those values are used rather than his
% recommended cp_liq - cp_vap = 2180, which yields cp_vap = 2040.

% Regarding ice, Ambaum mentions that the derivation assumes "the heat
% capacity at constant pressure is constant over the temperature range of
% interest" which is not good for ice, but they start with that assumption
% then check the result and find good accuracy. 




% From Romps:
% Rv = 461;                      % specific gas constant for water vapor [J kg-1 K-1]
% cp_liq = 4119;                 % specific heat capacity of liquid water
% cp_ice = 1861;
% cv_vap = 1418;                 % specific heat capacity of water vapor at constant volume [J kg-1 K-1]
% cv_liq = cp_liq;               % specific heat capacity of water at constant volume [J kg-1 K-1]
% cv_ice = cp_ice;               % specific heat capacity of ice at constant volume [J kg-1 K-1]
% cp_vap = cv_vap + Rv;          % specific heat capacity of water vapor at constant pressure [J kg-1 K-1]
% Ptrip = 611.65;                % triple point vapor pressure

% Romps' expression for saturation vapor pressure:
% pv_ice = pv0 * (T / T0) .^ (cp_vap - cv_ice)

% Note, the specific heat capacity of water vapor at constant pressure is
% computed as:
% cp_vap = cv_vap + Rv
% where
% cv_vap = 1384.5
%
% % Constants
% R = 8.314; % J/mol.K, Universal gas constant
% molar_mass_water_vapor = 18.015;  % g/mol, Molar mass of water vapor
%
% % Molar heat capacity at constant volume for water vapor (non-linear
% triatomic molecule):
% cv_molar = 6 / 2 * R;
%
% % Convert to specific heat capacity (J/kg.K)
% cv_vap = cv_molar / (molar_mass_water_vapor / 1000); % molar mass to kg/mol

   


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:3437fefe]
%   data: {"dataType":"symbolic","outputData":{"name":"D_e_expr","value":"D_{\\textrm{e0}} \\,{{\\left(\\frac{T}{T_f }\\right)}}^n"}}
%---
%[output:4f3b72fd]
%   data: {"dataType":"symbolic","outputData":{"name":"e_s_expr","value":"e_{\\textrm{s0}} \\,{\\mathrm{e}}^{\\frac{b\\,{\\left(T-T_f \\right)}}{T+c}}"}}
%---
%[output:3c5d355f]
%   data: {"dataType":"symbolic","outputData":{"name":"rho_v_expr","value":"\\frac{e_{\\textrm{s0}} \\,{\\mathrm{e}}^{\\frac{b\\,{\\left(T-T_f \\right)}}{T+c}} }{R_v \\,T}"}}
%---
%[output:3d0af837]
%   data: {"dataType":"symbolic","outputData":{"name":"d_k_v_dT","value":"\\begin{array}{l}\n\\frac{D_{\\textrm{e0}} \\,L_v \\,n\\,{{\\left(\\frac{T}{T_f }\\right)}}^{n-1} \\,{\\left(\\rho_v \\,\\sigma_1 -\\frac{\\rho_v }{T}\\right)}}{T_f }-D_e \\,L_v \\,{\\left(\\rho_v \\,{\\left(\\frac{2\\,b}{{{\\left(T+c\\right)}}^2 }-\\frac{2\\,b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^3 }\\right)}-\\rho_v \\,{\\sigma_1 }^2 -\\frac{2\\,\\rho_v }{T^2 }+\\frac{2\\,\\rho_v \\,\\sigma_1 }{T}\\right)}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{b}{T+c}-\\frac{b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^2 }\n\\end{array}"}}
%---
%[output:8722b300]
%   data: {"dataType":"symbolic","outputData":{"name":"d_k_v_dT","value":"\\begin{array}{l}\n\\frac{D_e \\,L_v \\,n\\,{\\left(\\rho_v \\,\\sigma_1 -\\frac{\\rho_v }{T}\\right)}}{T}-D_e \\,L_v \\,{\\left(\\rho_v \\,{\\left(\\frac{2\\,b}{{{\\left(T+c\\right)}}^2 }-\\frac{2\\,b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^3 }\\right)}-\\rho_v \\,{\\sigma_1 }^2 -\\frac{2\\,\\rho_v }{T^2 }+\\frac{2\\,\\rho_v \\,\\sigma_1 }{T}\\right)}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{b}{T+c}-\\frac{b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^2 }\n\\end{array}"}}
%---
%[output:0fecc6e0]
%   data: {"dataType":"symbolic","outputData":{"name":"d_e_s_dT","value":"e_s \\,{\\left(\\frac{b}{T+c}-\\frac{b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^2 }\\right)}"}}
%---
%[output:28891588]
%   data: {"dataType":"symbolic","outputData":{"name":"d_rho_v_dT","value":"\\rho_v \\,{\\left(\\frac{b}{T+c}-\\frac{b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^2 }\\right)}-\\frac{\\rho_v }{T}"}}
%---
%[output:2c83f7d8]
%   data: {"dataType":"symbolic","outputData":{"name":"d2_rho_v_dT2","value":"\\begin{array}{l}\n\\rho_v \\,{\\sigma_1 }^2 -\\rho_v \\,{\\left(\\frac{2\\,b}{{{\\left(T+c\\right)}}^2 }-\\frac{2\\,b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^3 }\\right)}+\\frac{2\\,\\rho_v }{T^2 }-\\frac{2\\,\\rho_v \\,\\sigma_1 }{T}\\\\\n\\mathrm{}\\\\\n\\textrm{where}\\\\\n\\mathrm{}\\\\\n\\;\\;\\sigma_1 =\\frac{b}{T+c}-\\frac{b\\,{\\left(T-T_f \\right)}}{{{\\left(T+c\\right)}}^2 }\n\\end{array}"}}
%---
%[output:0de36b50]
%   data: {"dataType":"symbolic","outputData":{"name":"e_s_expr","value":"T^c \\,a\\,{\\mathrm{e}}^{b\/T}"}}
%---
%[output:7e05bf15]
%   data: {"dataType":"symbolic","outputData":{"name":"de_s_dT_expr","value":"T^{c-1} \\,a\\,c\\,{\\mathrm{e}}^{b\/T} -\\frac{T^c \\,a\\,b\\,{\\mathrm{e}}^{b\/T} }{T^2 }"}}
%---
%[output:1bd7c419]
%   data: {"dataType":"symbolic","outputData":{"name":"de_s_dT","value":"\\frac{c\\,e_s }{T}-\\frac{b\\,e_s }{T^2 }"}}
%---
%[output:85d99bbb]
%   data: {"dataType":"symbolic","outputData":{"name":"d2e_s_dT2","value":"\\frac{b^2 \\,e_s }{T^4 }+\\frac{2\\,b\\,e_s }{T^3 }-\\frac{2\\,b\\,c\\,e_s }{T^3 }+\\frac{c\\,e_s \\,{\\left(c-1\\right)}}{T^2 }"}}
%---
%[output:7877f4b3]
%   data: {"dataType":"symbolic","outputData":{"name":"ans","value":"-\\frac{b\\,{\\mathrm{e}}^{b\/T} }{T^2 }"}}
%---
