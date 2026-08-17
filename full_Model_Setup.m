%% ========================================================================
%  PART I - FULL MODEL SETUP
%
%  This script prepares all data required for the numerical computation:
%    - furnace geometry and topology,
%    - angular quadrature,
%    - material properties,
%    - boundary conditions,
%    - initial temperature fields,
%    - slab and radiation contexts,
%    - inter-zone mass-flow and gas energy-balance models.
%
% ========================================================================

clear;
clc;


%% ========================================================================
%  1. MODEL INPUTS
% ========================================================================

% Angular quadrature
Nth = 4;
Nph = 8;

% Initial gas temperature
T_init = 1300;      % K

% Geometry files
meshFile   = 'input/gmsh_reheating_furnace.msh';
markerFile = 'input/nonParticipatingMedia.msh';

% Model configuration
solveWalls = false;

% Zone, slab, physical-constant, and wall definitions
run('input/subZonesDefinition.m');

slabsBbox = containers.Map( ...
    'KeyType', 'int32', ...
    'ValueType', 'any');

run('input/bramestoreSlab.m');
run('input/physicalConstant.m');
run('input/wallDataFurnace_Generation.m');


%% ========================================================================
%  2. GEOMETRY AND TOPOLOGY
% ========================================================================

mesh   = read_mesh(meshFile);
marker = read_marker_nonParticipatingMedia(markerFile);

[geom, mesh] = build_topology(mesh, marker);

slabs = build_slabs(slabsBbox);
clear slabsBbox

fprintf('Updating mesh and topology for slabs...\n');

[mesh, geom] = build_topology_slabs(mesh, geom, slabs);
mesh         = build_topology_zones(mesh, geom, zonesDefinition);

zonesCtx = build_zones_context( ...
    zonesDefinition, mesh, geom);


%% ========================================================================
%  3. ANGULAR QUADRATURE
% ========================================================================

ang = build_angles(Nth, Nph);

D = precompute_face_direction_weights( ...
    geom, ang);


%% ========================================================================
%  4. MATERIAL PROPERTIES
% ========================================================================

% Slabs
matProps.slabs.k = [ ...
     20.6397 ...
    -0.0414559 ...
     0.000103482 ...
    -7.14584e-8 ...
     1.59468e-11];

matProps.slabs.cp = [ ...
     256.912 ...
     0.885133 ...
     0.00306974 ...
    -1.98127e-5 ...
     3.87382e-8 ...
    -3.32115e-11 ...
     1.09459e-14];

matProps.slabs.rho = 7900;
matProps.slabs.eps = 1;


% Participating gas
matProps.fluid.kappa = 0.21;


% External-wall materials
matProps.extWall.material_ram_20_25.k = ...
    [2.58784 -0.00318676 1.80e-06];

matProps.extWall.material_ram_20_25.cp  = 500;
matProps.extWall.material_ram_20_25.rho = 2400;


matProps.extWall.material_flumisol.k = ...
    [0.0811675 -0.00022528 4.80e-07];

matProps.extWall.material_flumisol.cp  = 500;
matProps.extWall.material_flumisol.rho = 180;


matProps.extWall.material_g80f.k = ...
    [2.98951 0.00174679 -1.38e-06];

matProps.extWall.material_g80f.cp  = 500;
matProps.extWall.material_g80f.rho = 2770;


matProps.extWall.material_refor30.k = ...
    [0.487083 0.000658433 -3.16e-07];

matProps.extWall.material_refor30.cp  = 500;
matProps.extWall.material_refor30.rho = 2080;


matProps.extWall.material_eguiv_rl117_diatomee.k = ...
    [0.149595 -9.62e-06 0];

matProps.extWall.material_eguiv_rl117_diatomee.cp  = 500;
matProps.extWall.material_eguiv_rl117_diatomee.rho = 2400;


matProps.extWall.steel.k   = [16.27 0 0];
matProps.extWall.steel.cp  = 502.48;
matProps.extWall.steel.rho = 8030;


% Gas heat capacities
matProps.fluid.o2.cp = [ ...
     975.974 ...
    -0.66073 ...
     0.00223857 ...
    -2.451160e-06 ...
     1.281630e-09 ...
    -3.236540e-13 ...
     3.172820e-17];

matProps.fluid.n2.cp = [ ...
     1116.63 ...
    -0.615175 ...
     0.0015183 ...
    -1.296410e-06 ...
     5.498740e-10 ...
    -1.167650e-13 ...
     9.910130e-18];

matProps.fluid.h2o.cp = [ ...
     1875.84 ...
    -0.419047 ...
     0.00157042 ...
    -9.671970e-07 ...
     2.533150e-10 ...
    -2.497410e-14 ...
     0];

matProps.fluid.co2.cp = [ ...
     452.206 ...
     1.6783 ...
    -0.00140075 ...
     6.435410e-07 ...
    -1.534180e-10 ...
     1.475850e-14 ...
     0];


%% ========================================================================
%  5. BOUNDARY CONDITIONS
% ========================================================================

fprintf('Reading boundary conditions...\n');

% Slab conduction
bcConfig.slabs.symmetry = "adiabatic";
bcConfig.slabs.side     = "radiative";
bcConfig.slabs.top      = "radiative";
bcConfig.slabs.bottom   = "radiative";
bcConfig.slabs.front    = "radiative";
bcConfig.slabs.rear     = "dirichlet";

bcConfig.slabs.firstSlabRearT = 1533;     % K


% Radiation
bcConfig.rad.symmetry            = "symmetry";
bcConfig.rad.intWall             = "radiative";
bcConfig.rad.interfaceBlockedOff = "radiative";
bcConfig.rad.hs                  = "imposed_heat_flux";
bcConfig.rad.skid                = "imposed_heat_flux";

radBcType2Tags = int32([ ...
    886:900, ...
    1006:1053]);

radBcType3Tags = int32(946:1005);

q_imposed = 30500;     % W/m^2


%% ========================================================================
%  6. SOLVER PARAMETERS
% ========================================================================

% Coupled outer iterations
maxOuterIter = 30;
outerTol     = 1e-3;


% Slab conduction solver
solverParams = struct();

solverParams.slabs.picardTol     = 1e-4;
solverParams.slabs.picardMaxIter = 30;
solverParams.slabs.relaxAlpha    = 0.9;
solverParams.slabs.verbose       = true;


% Wall conduction solver
wallOpts = struct();

wallOpts.qConvMap = [];
wallOpts.qCorrMap = [];

wallOpts.solverParams = struct( ...
    'relaxAlpha', 0.9, ...
    'picardTol',  1e-6);

wallOpts.verbose = false;


% Gas energy-balance solver
zoneEnergyOpts = struct( ...
    'enforceMassConservation', true, ...
    'verbose', true);


%% ========================================================================
%  7. INITIAL TEMPERATURE FIELDS
% ========================================================================

radFields.T = T_init * ones( ...
    numel(mesh.cells.id), 1);

load('input/postprocZones_New1.mat');

zoneFields.T = postprocZones(:,2);

radFields = assign_zone_temperature_to_radFields( ...
    radFields, zoneFields, mesh, zonesCtx);


% Reference wall and skid data
load('input/extWallReport_New1.mat');
load('input/skid.mat');

Twall_zone_table = build_Twall_zone_table_from_reports( ...
    extWallReport, skid);


% Door temperatures imposed from the reference CFD data
Twall_zone_table('421_17') = 1500.7118;
Twall_zone_table('431_17') = 1481.1887;
Twall_zone_table('411_17') = 1483.6707;

Twall_zone_table('422_17') = 1500.7118;
Twall_zone_table('432_17') = 1481.1887;
Twall_zone_table('412_17') = 1483.6707;


% Inter-zone flow data
load('input/intSurfReport2.mat');


%% ========================================================================
%  8. SLAB CONDUCTION CONTEXT
% ========================================================================

fprintf('Preparing slabs for conduction...\n');

emptyCtx = struct( ...
    'id',             [], ...
    'Nc',             [], ...
    'globalCellIdx',  [], ...
    'globalToLocal',  [], ...
    'cellNeighbor',   [], ...
    'cellFace',       [], ...
    'faceBC',         [], ...
    'faceSubzone',    [], ...
    'cellVolume',     [], ...
    'cellCenter',     [], ...
    'faceArea',       [], ...
    'faceDist',       [], ...
    'surfaceFaces',   []);

slabCtxs = repmat( ...
    emptyCtx, numel(slabs), 1);

for i = 1:numel(slabs)

    slabCtxs(i) = build_slab_context( ...
        slabs(i), mesh, geom, bcConfig);

end

slabCtxs = build_slab_rear_coupling( ...
    slabCtxs, geom, bcConfig);


%% ========================================================================
%  9. RADIATION AND WALL CONTEXTS
% ========================================================================

radCtx = build_rad_context( ...
    geom, mesh);

radSurf = build_radiation_surfaces( ...
    zonesCtx, mesh, geom);


% Slab cells used for convergence monitoring
slabCells   = find(mesh.cells.isLoad);
T_prev_slab = radFields.T(slabCells);


% Wall conduction context
wallTargetTags = int32([ ...
    17 ...
    18 ...
    16 ...
    19 ...
    865 ...
    866]);

wallCtx = build_wall_context( ...
    wallDataFurnace, ...
    zonesCtx, ...
    geom, ...
    matProps, ...
    wallTargetTags);

wallState = [];
Tsi_prev  = [];


%% ========================================================================
%  10. GAS ENERGY-BALANCE MODEL
%
%  The inter-zone mass flows are extracted from the CFD solution and used
%  to construct the zonal energy-balance network. No gas-to-wall convective
%  loss is imposed here.
% ========================================================================

flowData = intSurfArray(:, [1 2]);

mergeGroups = {};

zEnCtx = build_zone_energy_context( ...
    zonesCtx, ...
    radSurf, ...
    geom, ...
    flowData, ...
    [], ...
    matProps.fluid.kappa, ...
    struct('mergeGroups', {mergeGroups}));

zEnCtx.Qconv_ez(:) = 0;


% Gas composition: [O2 N2 H2O CO2]
gasComposition = [ ...
    0.04 ...
    0.74 ...
    0.144 ...
    0.076];

gasProps = build_gas_properties( ...
    gasComposition);


% Burner heat-release rates [W]
% Zones: [111 211 311 312 322 332]
Qburner = [ ...
    1007555.1 ...
    1007555.1 ...
    1007555.1 ...
    118823.75 ...
    237647.496 ...
    237647.496];


% Burner mass-flow rates [kg/s]
% Lateral burners: [111 211 311]
mdotBurnerLat = [ ...
    0.4272 ...
    0.4272 ...
    0.4272];

% Roof burners: [312 322 332]
mdotBurnerRoof = [ ...
    0.07396 ...
    0.14792 ...
    0.14792];


src = struct( ...
    'Qburner',        Qburner, ...
    'mdotBurnerLat',  mdotBurnerLat, ...
    'mdotBurnerRoof', mdotBurnerRoof);


% Boundary energy and mass flows
% Qinlet zones:   [121 131]
% mdotOut zones:  [111 112 122 132]
% mdotInlet zones:[121 131]

src.Qinlet = [ ...
    1568200.7 ...
    5103877.96];

src.mdotOut = [ ...
    1.4705 ...
    1.2192 ...
    1.5920 ...
    1.2603614];

src.mdotInlet = [ ...
    0.7896 ...
    2.9007];


%% ========================================================================
%  11. INITIAL GAS ENERGY-ZONE TEMPERATURES
%
%  The initial temperature of each energy zone is obtained from the current
%  temperature field of the gas cells belonging to that zone.
% ========================================================================

T_ez = zeros(zEnCtx.nEZ,1);
Vsum = zeros(zEnCtx.nEZ,1);

for p = 1:zEnCtx.nSub

    ge = zEnCtx.gasElemOfSub(p);
    c  = double(radSurf.elem(ge).cells);
    e  = zEnCtx.ezOfSub(p);

    T_ez(e) = T_ez(e) + sum(radFields.T(c));
    Vsum(e) = Vsum(e) + numel(c);

end

T_ez = T_ez ./ max(Vsum,1);

Tzone_prev = T_ez;
