function Twall_zone_table = build_Twall_zone_table_from_reports(extWallReport, skid)

    Twall_zone_table = containers.Map( ...
        'KeyType', 'char', ...
        'ValueType', 'double');

    %% 1) External wall report
    Tmat = extWallReport.staticTemperature;

    [zoneIds, tags, values] = find(Tmat);

    for k = 1:numel(values)

        zoneTag = zoneIds(k);
        faceTag = tags(k);
        Tw      = values(k);

        if Tw == 0
            continue
        end

        key = sprintf('%d_%d', zoneTag, faceTag);
        Twall_zone_table(key) = Tw;
    end

    %% 2) Skid support report
    refractoryToTags = containers.Map( ...
        'KeyType', 'double', ...
        'ValueType', 'any');

    refractoryToTags(10) = [925 926];
    refractoryToTags(1)  = 931;
    refractoryToTags(8)  = [932 936];
    refractoryToTags(9)  = 937;
    refractoryToTags(11) = [941 942];

    for iz = 1:numel(skid.zoneID)

        zoneTag = skid.zoneID(iz);

        for ir = 1:numel(skid.refractoryID)

            refractoryID = skid.refractoryID(ir);
            Tw = skid.static_temperature_K(iz, ir);

            if isnan(Tw) || Tw == 0
                continue
            end

            if ~isKey(refractoryToTags, refractoryID)
                continue
            end

            faceTags = refractoryToTags(refractoryID);

            for jt = 1:numel(faceTags)

                faceTag = faceTags(jt);
                key = sprintf('%d_%d', zoneTag, faceTag);

                Twall_zone_table(key) = Tw;
            end
        end
    end
end