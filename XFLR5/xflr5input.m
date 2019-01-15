filename = 'mainwingdata.xlsx'; % xflr5 data
sheet = 1;

XFLR5 = {};

% freestream air speed (m/s)
Vinf_set = 'B5';
XFLR5.Vinf = xlsread(filename,sheet,Vinf_set);

% angle of attack (deg)
alpha_set = 'B6';
XFLR5.alpha = xlsread(filename,sheet,alpha_set);

% lift coefficient of the wing
CL_wing_set = 'B10';
XFLR5.CL_wing = xlsread(filename,sheet,CL_wing_set);

% drag coefficient of the wing
CD_wing_set = 'B12';
XFLR5.CD_wing = xlsread(filename,sheet,CD_wing_set);

% bending coefficient of the wing
Cm_wing_set = 'B14';
XFLR5.Cm_wing = xlsread(filename,sheet,Cm_wing_set);

% x-position of the wing center of pressure
XCP_set = 'B16';
XFLR5.XCP = xlsread(filename,sheet,XCP_set);

% z-position of the wing center of pressure
ZCP_set = 'F16';
XFLR5.ZCP = xlsread(filename,sheet,ZCP_set);

% local span position (m)
span_position_set = 'A22:A61';
XFLR5.span_position = xlsread(filename,sheet,span_position_set);

% local chord length (m)
chord_set = 'B22:B61';
XFLR5.chord = xlsread(filename,sheet,chord_set);

% local induced angle of attack
alpha_i_set = 'C22:C61';
XFLR5.alpha_i = xlsread(filename,sheet,alpha_i_set);

% local lift coefficient
cl_set = 'D22:D61';
XFLR5.cl = xlsread(filename,sheet,cl_set);

% local profile drag coefficient
cd_p_set = 'E22:E61';
XFLR5.cd_p = xlsread(filename,sheet,cd_p_set);

% local induced drag coefficient
cd_i_set = 'F22:F61';
XFLR5.cd_i = xlsread(filename,sheet,cd_i_set);

% local pitching moment due to pressure and viscous forces
cm_geo_set = 'G22:G61';
XFLR5.cm_geo = xlsread(filename,sheet,cm_geo_set);

% local position center of pressure
xcp_set = 'K22:K61';
XFLR5.xcp = xlsread(filename,sheet,xcp_set);

% local bending moment
M_bend_set = 'L22:L61';
XFLR5.M_bend = xlsread(filename,sheet,M_bend_set);

% local transition laminar to turbulent flow
xtr_top_set = 'I22:I61';
XFLR5.xtr_top = xlsread(filename,sheet,xtr_top_set);
xtr_bottom_set = 'I22:I61';
XFLR5.xtr_bottom = xlsread(filename,sheet,xtr_bottom_set);

clearvars -except XFLR5
