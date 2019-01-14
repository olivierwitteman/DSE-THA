clear all

% load XFLR5 data
xflr5input

% modify these manually
S = 8.2861;                             % wing area (m^2)
span = 9.1;                             % wing span (m)
rho = 1.225;                            % air density (kg/m^3)
MAC = 0.996;                            % mean aerodynamic chord (m)
c_r = 1.3;                              % root chord (m)
c_t = 0.52;                             % tip chord (m)

stepi = 0.01;                           % spacing for interpolation (m)

% from the data sheet
Vinf = XFLR5.Vinf;                      % freestream air speed (m/s)
q = 0.5*rho*Vinf^2;                     % dynamic pressure (Pa)
span_position = XFLR5.span_position;    % span position (m)
chord = XFLR5.chord;                    % local chord length

% grid vector along the span (m)
gridf = (-span/2:stepi:span/2)';
gridi = gridf(gridf>=min(span_position) & gridf<=max(span_position));

% local chord along grid vector (m)
chordi = interp1(span_position,chord,gridi);

% local lift per unit span (N/m)
cl = XFLR5.cl;
L = interp1(span_position,q*cl.*chord,gridi);
L = [0; L; 0];

% local profile drag per unit span (N/m)
cd_p = XFLR5.cd_p;
D_p = interp1(span_position,q*cd_p.*chord,gridi);
D_p_first = interp1(gridi,D_p,gridf(1),'linear','extrap');
D_p_end = interp1(gridi,D_p,gridf(end),'linear','extrap');
D_p = [D_p_first; D_p; D_p_end];

% local induced drag per unit span (N/m)
cd_i = XFLR5.cd_i;
D_i = interp1(span_position,q*cd_i.*chord,gridi);
D_i_first = interp1(gridi,D_i,gridf(1),'linear','extrap');
D_i_end = interp1(gridi,D_i,gridf(end),'linear','extrap');
D_i = [D_i_first; D_i; D_i_end];

% local drag per unit span (N/m)
cd = cd_i + cd_p;
D = D_i + D_p;

% local pitching moment due to pressure and viscous forces (N*m)
cm_geo = XFLR5.cm_geo;
M_pitch = interp1(span_position,q*S*chord.*cm_geo,gridi);
M_pitchfirst = interp1(gridi,M_pitch,gridf(1),'linear','extrap');
M_pitchend = interp1(gridi,M_pitch,gridf(end),'linear','extrap');
M_pitch = [M_pitchfirst; M_pitch; M_pitchend];
clear M_pitchfirst M_pitchend

% local bending moment (N*m)
M_bend = XFLR5.M_bend;
M_bend_i = interp1(span_position,M_bend,gridi);
M_bend_i = [0; M_bend_i; 0];

% local position center of pressure
xcp = XFLR5.xcp;