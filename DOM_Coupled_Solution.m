%% ========================================================================
%  PART II - FULL DOM COUPLED SOLUTION
%
%  This script solves the coupled radiation / wall conduction / slab
%  conduction problem using the full Discrete Ordinates Method (DOM).
%
% ========================================================================

I = [];


for outerIter = 1:maxOuterIter

    %% Radiation properties


    radProps = build_radiation_properties( ...
        geom, ...
        mesh, ...
        radFields, ...
        Twall_zone_table, ...
        matProps.fluid.kappa, ...
        physicalConst, ...
        radBcType2Tags, ...
        radBcType3Tags);



    %% DOM radiation solution

    solverParams.I0 = I;


    [I, infoRad, radProps] = solve_radiation( ...
        geom, ...
        mesh, ...
        radProps, ...
        radCtx, ...
        D, ...
        ang, ...
        q_imposed, ...
        physicalConst, ...
        solverParams);



    %% Incident radiation field


    [G_face, ~] = compute_G_field( ...
        I, ...
        geom, ...
        mesh, ...
        radProps, ...
        D, ...
        ang);

    radFields.G = G_face;



    %% Wall conduction

    if solveWalls

        wallOpts.wallState = wallState;

        [Twall_zone_table, wallState, infoWall] = ...
            solve_walls_conduction( ...
                wallCtx, ...
                radFields, ...
                geom, ...
                Twall_zone_table, ...
                physicalConst, ...
                wallOpts);

    end


    %% Slab conduction

    [radFields, slabCtxs, infoCond] = ...
        solve_slabs_conduction( ...
            slabCtxs, ...
            radFields, ...
            matProps, ...
            "steady", ...
            [], ...
            solverParams);


    %% Coupled convergence

    T_new_slab = radFields.T(slabCells);

    dT = norm(T_new_slab - T_prev_slab) / ...
         max(norm(T_new_slab), eps);


    if solveWalls

        Tsi_new = infoWall.Tsi;

        if isempty(Tsi_prev)

            dTwall = inf;

        else

            dTwall = norm(Tsi_new - Tsi_prev) / ...
                     max(norm(Tsi_new), eps);

        end

        wallImbal = infoWall.maxAbsImbalance;

    else

        Tsi_new   = [];
        dTwall    = 0;
        wallImbal = 0;

    end


    if dT < outerTol && ...
       dTwall < outerTol && ...
       outerIter > 1

        fprintf('  -> Coupled solution converged.\n');
        break

    end


    T_prev_slab = T_new_slab;
    Tsi_prev    = Tsi_new;