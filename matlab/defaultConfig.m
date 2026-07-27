function params = defaultConfig()
    params = struct();

    % Experiment metadata
    params.experimentId = "baseline_sdpt3";
    params.caseId = "case_fixed_001";
    params.algorithm = "rss";
    params.solver = "sdpt3";

    % Robot geometry
    params.Lx = 0.655;
    params.Ly = 0.335;
    params.a = params.Lx / 2;
    params.b = params.Ly / 2;
    params.wheel_pos = [
        params.a,  params.b;
        -params.a, params.b;
        -params.a, -params.b;
        params.a,  -params.b
    ];

    params.vimax = 5;
    params.phidotmax = 5 * pi;

    % Legacy controller parameters
    params.k1 = 0.15;
    params.k2 = 0.15;
    params.k3 = 0.1;
    params.eps = 0.001;

    % Simulation
    params.dt = 0.01;
    params.t_end = 1.0;
    params.num_steps = round(params.t_end / params.dt);

    % Reference
    params.ctrl_pts = [
        0,     0;
        0.750, 0.250;
        0.250, 0.750;
        1.250, -1.000;
        1.000, 0
    ];
    params.num_path_pts = params.num_steps;

    % Protocol
    params.trajectory.mode = "paper_fixed";
    params.trajectory.timeProtocol = "legacy_100_points";

    % Solver selection
    params.use_hpipm = false;

    % Output
    params.output.livePlot = false;
    params.output.saveFigures = true;
    params.output.saveFullLog = true;
    params.output.verbose = true;
end
