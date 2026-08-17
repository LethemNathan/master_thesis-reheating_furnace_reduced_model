
R = 1.485; %m

L = 6.490;
H = 3.7;
x0 = -29.2095;
x1 = -37.275;
slabHeight = 2.3;

zonesDefinition = struct();

zonesDefinition.R = R;



zonesDefinition.volumes  = struct(...
        'zoneTag',{},...
        'bbox',{});

volumes  = struct(...
        'zoneTag',{},...
        'bbox',{});

volumes(end+1).zoneTag = int32(111);
volumes(end).bbox = [-31.669 L-R 0 x0 L slabHeight];

volumes(end+1).zoneTag = int32(112);
volumes(end).bbox = [-31.669 L-R slabHeight x0 L H];

volumes(end+1).zoneTag = int32(121);
volumes(end).bbox = [-31.669 (L-R)/2 0 x0 L-R slabHeight];

volumes(end+1).zoneTag = int32(122);
volumes(end).bbox = [-31.669 (L-R)/2 slabHeight x0 L-R H];

volumes(end+1).zoneTag = int32(131);
volumes(end).bbox = [-31.669 0 0 x0 (L-R)/2 slabHeight];

volumes(end+1).zoneTag = int32(132);
volumes(end).bbox = [-31.669 0 slabHeight x0 (L-R)/2 H];




volumes(end+1).zoneTag = int32(211);
volumes(end).bbox = [-34.140 L-R 0 -31.669 L slabHeight];

volumes(end+1).zoneTag = int32(212);
volumes(end).bbox = [-34.140 L-R slabHeight -31.669 L H];

volumes(end+1).zoneTag = int32(221);
volumes(end).bbox = [-34.140 (L-R)/2 0 -31.669 L-R slabHeight];

volumes(end+1).zoneTag = int32(222);
volumes(end).bbox = [-34.140 (L-R)/2 slabHeight -31.669 L-R H];

volumes(end+1).zoneTag = int32(231);
volumes(end).bbox = [-34.140 0 0 -31.669 (L-R)/2 slabHeight];

volumes(end+1).zoneTag = int32(232);
volumes(end).bbox = [-34.140 0 slabHeight -31.669 (L-R)/2 H];




volumes(end+1).zoneTag = int32(311);
volumes(end).bbox = [-36.132 L-R 0 -34.140 L slabHeight];

volumes(end+1).zoneTag = int32(312);
volumes(end).bbox = [-36.132 L-R slabHeight -34.140 L H];

volumes(end+1).zoneTag = int32(321);
volumes(end).bbox = [-36.132 (L-R)/2 0 -34.140 L-R slabHeight];

volumes(end+1).zoneTag = int32(322);
volumes(end).bbox = [-36.132 (L-R)/2 slabHeight -34.140 L-R H];

volumes(end+1).zoneTag = int32(331);
volumes(end).bbox = [-36.132 0 0 -34.140 (L-R)/2 slabHeight];

volumes(end+1).zoneTag = int32(332);
volumes(end).bbox = [-36.132 0 slabHeight -34.140 (L-R)/2 H];




volumes(end+1).zoneTag = int32(411);
volumes(end).bbox = [x1 L-R 0 -36.132 L slabHeight];

volumes(end+1).zoneTag = int32(412);
volumes(end).bbox = [x1 L-R slabHeight -36.132 L H];

volumes(end+1).zoneTag = int32(421);
volumes(end).bbox = [x1 (L-R)/2 0 -36.132 L-R slabHeight];

volumes(end+1).zoneTag = int32(422);
volumes(end).bbox = [x1 (L-R)/2 slabHeight -36.132 L-R H];

volumes(end+1).zoneTag = int32(431);
volumes(end).bbox = [x1 0 0 -36.132 (L-R)/2 slabHeight];

volumes(end+1).zoneTag = int32(432);
volumes(end).bbox = [x1 0 slabHeight -36.132 (L-R)/2 H];


zonesDefinition.volumes = volumes;

zonesDefinition.numZones = numel(zonesDefinition.volumes);

zonesDefinition.lateralBurnerLoc = [111 211 311];
zonesDefinition.roofBurnerLoc = [312 322 332];
%zonesDefinition.inletLoc = [111 112 121 122 131 132 711 712 721 722 732 731];
zonesDefinition.inletLoc = [111 112 121 122 131 132];

clear volumes
clear L
clear H
clear x0
clear x1
clear R
clear slabHeight





