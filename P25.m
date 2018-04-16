clc
clear

format short g

%% install toolboxes
% only run if the toolbox is not installed
if isempty(ver('GEOM3D'))
    toolboxFile = 'geom3d-2017.12.01.mltbx';
    installedToolbox = matlab.addons.toolbox.installToolbox(toolboxFile);
end

if isempty(ver('GEOM2D'))
    toolboxFile = 'geom2d-2017.08.31.mltbx';
    installedToolbox = matlab.addons.toolbox.installToolbox(toolboxFile);
end

%% parmeters
p.g     = 10;
p.rho_w = 1;
p.rho_b = 0.4*p.rho_w;

p.l = 3;
p.w = 5; 
p.h = 1;

p.m = p.l*p.w*p.h*p.rho_b;

p.I_B = p.m/12*[(p.l^2 + p.h^2) 0               0;
                 0              (p.w^2 + p.h^2) 0;
                 0              0               (p.l^2 + p.w^2)];

%% Part A equilibrium point
p.hw = (p.rho_b / p.rho_w) * p.h;

%% Part Bi - Euler's Bouyancy Method
sub_vol = p.l*p.w*p.hw;
theta  = deg2rad(10);
% this is a triangle
M_W_B  = [1/12*p.rho_w*p.g*p.l*theta*(p.w)^3; 0; 0];
r_C0relW_B = [0; 0; -p.hw/2];
r_GrelW_B  = [0; 0;  p.h/2-p.hw];
F_grav_B   = p.m*p.g*[0; sin(theta); -cos(theta)];
F_bouy_B   = sub_vol*p.rho_w*p.g*[0; -sin(theta); cos(theta)];
Msum_W     = M_W_B + cross(r_C0relW_B, F_bouy_B) + cross(r_GrelW_B, F_grav_B);

%% initial conditions
boat_G_F_0  = [0; 0; p.h/2 - 0.8*p.hw];

v_G_F_0     = [0; 0; 0];

R_0 = angle2dcm(deg2rad(0), deg2rad(-4), theta);

w_B_0 = [0; 0; 0];  w_F_0 = R_0 * w_B_0;

z_F_0 = [boat_G_F_0; v_G_F_0; w_F_0; reshape(R_0, 9,1)];

%% create boat and water
% boat center of gravity and rotation initialization

lx = p.w;
ly = p.l;
lz = p.h;
p.boat_v_relG_B = [ 0  0  0;  0  0 lz;  0 ly lz;  0 ly  0;  % 1-4
                    0 ly  0;  0 ly lz; lx ly lz; lx ly  0;  % 5-8
                   lx  0  0; lx  0 lz; lx ly lz; lx ly  0;  % 9-12
                    0  0  0;  0  0 lz; lx  0 lz; lx  0  0;  % 13-16
                    0  0 lz;  0 ly lz; lx ly lz; lx  0 lz;  % 17-20
                    0  0  0;  0 ly  0; lx ly  0; lx  0  0] - [0.5*lx 0.5*ly 0.5*lz];     
f           = reshape(1:24, 4,6)';
v_all_F_0   = trans_pts(rot_pts(p.boat_v_relG_B, R_0), boat_G_F_0);

figure(1)
clf
view([52, 33])
grid on
hold on

% draw boat
boat = patch('Faces', f, 'Vertices', v_all_F_0, 'FaceColor', 'r');  
set(boat, 'facealpha', 0.4);

% draw center of mass
COM  = plot3(boat_G_F_0(1), boat_G_F_0(2), boat_G_F_0(3), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 8);

% find submerged part
sub_verts = findSubmerged(boat);
waterline = plot3(sub_verts(:,1), sub_verts(:,2), sub_verts(:,3), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% point of application of submerged body
sub_C_F = polyhedronCentroid(sub_verts, boat.Faces);
sub_C   = plot3(sub_C_F(1), sub_C_F(2), sub_C_F(3), 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

% draw water
xWat  = [-100   100  100  -100  ];
yWat  = [-100  -100  100   100  ];
zWat  = [ p.hw  p.hw p.hw  p.hw ];
water       = patch(xWat, yWat, zWat, 'FaceColor', 'b');
set(water,'facealpha',0.2)
set(water,'edgealpha',0.2)

% more graphing stuff
axis 'equal'
offset = [boat_G_F_0(1) boat_G_F_0(1) boat_G_F_0(2) boat_G_F_0(2) boat_G_F_0(3) boat_G_F_0(3)];
axis(offset + [-lx lx -ly ly -2*lz lz])
xlabel('X Axis');   ylabel('Y Axis');   zlabel('Z Axis');

%% run equation solver
tf  = 20;
n   = tf*120;
t   = linspace(0, tf, n);

opts.RelTol = 1E-5;
opts.AbsTol = 1E-5;
[t, sol] = ode45(@RHS_cent, t, z_F_0, opts, p);

r_G_sol = sol(:, 1:3);
v_G_sol = sol(:, 4:6);
w_sol   = sol(:, 7:9);
R_sol   = sol(:, 10:18);

%% animate results
filename = 'Buoyant Block';
v = VideoWriter(filename);
open(v);
time_frac = 1;
tic
start = toc;
i = 0;
while start < tf
    i = i+1;
    % timing stuff
    start = toc*time_frac;
    [d, ind_now] = min(abs(t-start));
    
    
    % grab new info
    boat_G_F     = r_G_sol(ind_now,:);
    R            = reshape(R_sol(ind_now, :), 3,3);
    w_F          = w_sol(ind_now, :);
    
    % update boat points
    v_boat_F      = trans_pts(rot_pts(p.boat_v_relG_B, R), boat_G_F');
    boat.Vertices = v_boat_F;
    
    % update COM
    COM.XData = boat_G_F(1);
    COM.YData = boat_G_F(2);
    COM.ZData = boat_G_F(3);
    
    % update submerged
    sub_verts       = findSubmerged(boat);
    waterline.XData = sub_verts(:,1);
    waterline.YData = sub_verts(:,2);
    waterline.ZData = sub_verts(:,3);
    sub_C_F         = polyhedronCentroid(sub_verts, boat.Faces);
    sub_C.XData     = sub_C_F(1);
    sub_C.YData     = sub_C_F(2);
    sub_C.ZData     = sub_C_F(3);
    
    % recording stuff
%     writeVideo(v,getframe(gcf));
    
    % timing stuff
    dt = tf/n - (toc*time_frac - start);
    pause(dt);
end

%% functions
function zdot = RHS_cent(t, z_t, p)
    % ------------------ STATES -----------------------
    r_G_F = z_t(1:3);
    v_G_F = z_t(4:6);
    w_F   = z_t(7:9);
    R     = reshape(z_t(10:18), 3,3);
    
    
    % ---------------- GEOMETRY -----------------------
    boat.Vertices   = trans_pts(rot_pts(p.boat_v_relG_B, R), r_G_F);
    f               = reshape(1:24, 4,6)';
    sub_verts       = sortrows(findSubmerged(boat));
    [K, vol_sub]    = convhull(sub_verts(:,1), sub_verts(:,2), sub_verts(:,3));
    sub_C_F         = polyhedronCentroid(sub_verts, f)';
    
    % ---------------- DYNAMICS -----------------------
    % LMB    
    F_grav_F    = [0; 0; -p.g*p.m];
    F_bouy_F    = [0; 0;  p.g*p.rho_w*vol_sub];
    Fsum_F      = F_grav_F + F_bouy_F;
    a_G_F       = Fsum_F/p.m;
    
    % AMB
    Msum_G_F    = cross(sub_C_F, F_bouy_F) + cross(r_G_F, F_grav_F);
    I           = R*p.I_B*R';
    w_dot       = I \ ( Msum_G_F - cross(w_F, I*w_F) );
    R_dot       = skew(w_F)*R;
    
    zdot        = [v_G_F; a_G_F; w_dot; reshape(R_dot, 9,1)];
end

% find submerged vertices
function sub_verts = findSubmerged(boat)
    for i = 1:4
        sub_verts(2*i-1,1:3) = boat.Vertices(20+i, 1:3);
        [sub_verts(2*i,1:3),check] = plane_line_intersect([0; 0; 1],zeros(1,3),boat.Vertices(16+i, 1:3),boat.Vertices(20+i, 1:3));
%         if check == 3
%             warning('Top or bottom of the boat has crossed water line for at least one corner. Expect erroneous results.');
%         end
    end
end

function pts_p = rot_pts(pts, R)
    pts_p = zeros(size(pts));
    for i = 1:size(pts,1)
        pts_p(i,1:3) = (R*pts(i,1:3)')';
    end
end

function pts_p = trans_pts(pts, d)
    pts_p = pts + d';
end
