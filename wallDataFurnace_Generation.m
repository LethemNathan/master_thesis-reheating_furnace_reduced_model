


%% Wall data for furnace thermal model

% Default external boundary conditions
hAmb_default      = 10;     % W/(m^2 K) - to be adapted
Tinf_default      = 293; % K - ambient/free stream temperature
TextRad_default   = 293; % K - external radiation temperature

% Default mesh density
nodesPerMeter_default = 100; % nodes/m - to be adapted

wallDataFurnace = struct([]);

%% Tag 18 - Right wall
wallDataFurnace(end+1).tag = 18;
wallDataFurnace(end).wallName = 'Right wall';
wallDataFurnace(end).nLayers = 2;
wallDataFurnace(end).hAmb = hAmb_default;
wallDataFurnace(end).Tinf = Tinf_default;
wallDataFurnace(end).TextRad = TextRad_default;

wallDataFurnace(end).layers(1).materialName = 'material_ram_20-25';
wallDataFurnace(end).layers(1).thickness = 0.255;
wallDataFurnace(end).layers(1).nodesPerMeter = nodesPerMeter_default;

wallDataFurnace(end).layers(2).materialName = 'material_flumisol';
wallDataFurnace(end).layers(2).thickness = 0.165;
wallDataFurnace(end).layers(2).nodesPerMeter = nodesPerMeter_default;


%% Tag 16 - Sole
wallDataFurnace(end+1).tag = 16;
wallDataFurnace(end).wallName = 'Sole';
wallDataFurnace(end).nLayers = 3;
wallDataFurnace(end).hAmb = 2.13;
wallDataFurnace(end).Tinf = 399;
wallDataFurnace(end).TextRad = TextRad_default;

wallDataFurnace(end).layers(1).materialName = 'material_g80f';
wallDataFurnace(end).layers(1).thickness = 0.116;
wallDataFurnace(end).layers(1).nodesPerMeter = nodesPerMeter_default;

wallDataFurnace(end).layers(2).materialName = 'material_refor30';
wallDataFurnace(end).layers(2).thickness = 0.132;
wallDataFurnace(end).layers(2).nodesPerMeter = nodesPerMeter_default;

wallDataFurnace(end).layers(3).materialName = 'material_eguiv_rl117_diatomee';
wallDataFurnace(end).layers(3).thickness = 0.262;
wallDataFurnace(end).layers(3).nodesPerMeter = nodesPerMeter_default;


%% Tag 19 - Roof
wallDataFurnace(end+1).tag = 19;
wallDataFurnace(end).wallName = 'Roof';
wallDataFurnace(end).nLayers = 2;
wallDataFurnace(end).hAmb = hAmb_default;
wallDataFurnace(end).Tinf = Tinf_default;
wallDataFurnace(end).TextRad = TextRad_default;

wallDataFurnace(end).layers(1).materialName = 'material_ram_20-25';
wallDataFurnace(end).layers(1).thickness = 0.070;
wallDataFurnace(end).layers(1).nodesPerMeter = nodesPerMeter_default;

wallDataFurnace(end).layers(2).materialName = 'material_flumisol';
wallDataFurnace(end).layers(2).thickness = 0.230;
wallDataFurnace(end).layers(2).nodesPerMeter = nodesPerMeter_default;


%% Tag 17 - Door
wallDataFurnace(end+1).tag = 17;
wallDataFurnace(end).wallName = 'Door';
wallDataFurnace(end).nLayers = 1;
wallDataFurnace(end).hAmb = hAmb_default;
wallDataFurnace(end).Tinf = Tinf_default;
wallDataFurnace(end).TextRad = TextRad_default;

wallDataFurnace(end).layers(1).materialName = 'steel';
wallDataFurnace(end).layers(1).thickness = 0.23; % unknown in the table
wallDataFurnace(end).layers(1).nodesPerMeter = nodesPerMeter_default;


%% Tag 865 - Roof
wallDataFurnace(end+1).tag = 865;
wallDataFurnace(end).wallName = 'Roof';
wallDataFurnace(end).nLayers = 2;
%wallDataFurnace(end).hAmb = 30;
%wallDataFurnace(end).Tinf = 358;
wallDataFurnace(end).hAmb = hAmb_default;
wallDataFurnace(end).Tinf = Tinf_default;
wallDataFurnace(end).TextRad = TextRad_default;

wallDataFurnace(end).layers(1).materialName = 'material_ram_20-25';
wallDataFurnace(end).layers(1).thickness = 0.070;
wallDataFurnace(end).layers(1).nodesPerMeter = nodesPerMeter_default;

wallDataFurnace(end).layers(2).materialName = 'material_flumisol';
wallDataFurnace(end).layers(2).thickness = 0.230;
wallDataFurnace(end).layers(2).nodesPerMeter = nodesPerMeter_default;


%% Tag 866 - Roof
wallDataFurnace(end+1).tag = 866;
wallDataFurnace(end).wallName = 'Roof';
wallDataFurnace(end).nLayers = 2;
wallDataFurnace(end).hAmb = hAmb_default;
wallDataFurnace(end).Tinf = Tinf_default;
wallDataFurnace(end).TextRad = TextRad_default;

wallDataFurnace(end).layers(1).materialName = 'material_ram_20-25';
wallDataFurnace(end).layers(1).thickness = 0.070;
wallDataFurnace(end).layers(1).nodesPerMeter = nodesPerMeter_default;

wallDataFurnace(end).layers(2).materialName = 'material_flumisol';
wallDataFurnace(end).layers(2).thickness = 0.230;
wallDataFurnace(end).layers(2).nodesPerMeter = nodesPerMeter_default;