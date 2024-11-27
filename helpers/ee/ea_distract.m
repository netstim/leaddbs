function ea_distract()
% Lead-DBS Neural Navigator: A minimal Easter Egg game within Lead-DBS controlled by the space bar.
% The player navigates a DBS electrode through neural pathways, avoiding obstacles.

%% Global Game Variables
space_pressed = false;   % Space bar press flag
zaps_remaining = 3;       % Player starts with 3 zaps
zap_pressed = false;     % Zap button press flag
isGameOver = false;      % Game over flag
waitForRestart = false;  % Restart flag
fig = [];                % Figure handle
mscCredits = [];         % Music Credit Text
if exist('splat.mat', 'file')
    splatSound = load('splat.mat');
else
    warning('splat.mat not found. Sound effects will be disabled.');
    splatSound = [];
end

if exist('handel.mat', 'file')
    handelSound = load('handel.mat');
else
    warning('handel.mat not found. Sound effects will be disabled.');
    handelSound = [];
end

if exist('gong.mat', 'file')
    gongSound = load('gong.mat');
else
    warning('handel.mat not found. Sound effects will be disabled.');
    gongSound = [];
end

%% Set Up Figure Window
fig = figure('Name', 'Lead-DBS Neural Navigator', 'NumberTitle', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none', 'Color', 'w', ...
    'KeyPressFcn', @keyDown, 'KeyReleaseFcn', @keyUp, ...
    'CloseRequestFcn', @closeGame, ...
    'WindowState', 'maximized');  % Maximize the window on launch

pause(0.5);

% Set the axes to fill the figure window, leaving space at the top for the title
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0, 0, 1, 0.95]);

% Hide axes ticks and labels but keep the title visible
set(ax, 'XTick', [], 'YTick', [], 'Box', 'off');
hold(ax, 'on');

% Update axes limits to fill the screen while maintaining aspect ratio
updateAxesLimits();

% Handle window resizing
set(fig, 'ResizeFcn', @(src, event) updateAxesLimits());


%% Create and Play Background Music
% Create the melody
bgsound=load(fullfile(ea_getearoot,'helpers','ee','bg.mat'));

% Create an audioplayer object
player = audioplayer(bgsound.melody, bgsound.fs);

% Set up the audioplayer to loop the melody
set(player, 'StopFcn', @(~,~) play(player));

% Play the melody in the background
play(player);


%% Start Game Loop
while true
    % Display Intro Screen
    introScreen();

    % Wait for player to start
    waitForStart = true;
    set(fig, 'KeyPressFcn', @startKeyDown);
    while waitForStart
        pause(0.1);
        if ~ishandle(fig)
            return; % Exit if figure is closed
        end
    end

    % Start Countdown
    countdownSequence();

    % Reset Game Variables
    [score, isGameOver] = runGame();

    % Game Over Screen
    if ishandle(fig)
        ea_chirp; ea_chirp; ea_chirp;

        pause(3); % make sure player stops pressing buttons
        % Fetch commit data from GitHub
        commitMessages = fetchCommitData();

        % Set waitForRestart to true before displaying end credits
        waitForRestart = true;
        set(fig, 'KeyPressFcn', @restartKeyDown);

        cla(ax);
        xLim = get(ax, 'XLim');
        yLim = get(ax, 'YLim');
        centerX = mean(xLim);
        text(centerX, 60, 'Game Over!', ...
            'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Parent', ax);
        text(centerX, 50, ['Fiber Score: ', num2str(score)], ...
            'Color', 'k', 'FontSize', 16, 'HorizontalAlignment', 'center', 'Parent', ax);

        text(centerX, 20, 'Press Space to Restart or Esc to Exit', ...
            'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Parent', ax);



        % Display scrolling end credits
        displayEndCredits(commitMessages');

        % Wait for space bar to restart or escape to exit
        while waitForRestart
            showMusicCredits;
            if ~ishandle(fig)
                return; % Exit if figure is closed
            end
            pause(0.5)
            delMusicCredits;
        end
        zaps_remaining = 3; % reset remaining zaps
    else
        break; % Exit the game loop if figure is closed
    end
end

%% Nested Functions

    function introScreen()
        % Display the intro screen with the provided ASCII art
        cla(ax);
        % Get current axes limits
        xLim = get(ax, 'XLim');
        yLim = get(ax, 'YLim');
        centerX = mean(xLim);
        centerY = mean(yLim)+(yLim(2)-yLim(1))/10;

        ascii_art = {
            '                                                Welcome to                                               ', ...
            '                                                                                                         ', ...
            '                  _______________________________________                                                ', ...
            '                 /           ______     ___      _____   \                                               ', ...
            '                 | | ocali-  |         /   \    |     \  |             _______                           ', ...
            '                 | | zation  | lec     |   |    | BS   \ |          = / O   / \                          ', ...
            '                 | |         |---      /___\    | Sur- | |         = (   \./___)                         ', ...
            '                 | | of      | trodes |     |   | gery / |        =   \_ __ _)                           ', ...
            '                 | |_____    |_____   / fter \  |_____/  |         =   (O)(O)                            ', ...
            '                 \_______________________________________/                                               ', ...
            '                                                                                                         ', ...
            '                                                 * * * * *                                               ', ...
            '                                                                                                         ', ...
            ' #     #                                       #     #                                                   ', ...
            ' ##    # ###### #    # #####    ##   #         ##    #   ##   #    # #  ####    ##   #####  ####  #####  ', ...
            ' # #   # #      #    # #    #  #  #  #         # #   #  #  #  #    # # #    #  #  #    #   #    # #    # ', ...
            ' #  #  # #####  #    # #    # #    # #         #  #  # #    # #    # # #      #    #   #   #    # #    # ', ...
            ' #   # # #      #    # #####  ###### #         #   # # ###### #    # # #  ### ######   #   #    # #####  ', ...
            ' #    ## #      #    # #   #  #    # #         #    ## #    #  #  #  # #    # #    #   #   #    # #   #  ', ...
            ' #     # ######  ####  #    # #    # ######    #     # #    #   ##   #  ####  #    #   #    ####  #    # '
            };
        text(centerX, centerY, ascii_art, 'Color', 'k', 'FontName', 'FixedWidth', ...
            'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Parent', ax, 'Interpreter', 'none');

        text(centerX, (yLim(1) + (yLim(2)-yLim(1))/15), 'Press Space to Start and to control the Electrode - press . to zap!', ...
            'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Parent', ax, ...
            'Interpreter', 'none');
        drawnow;
    end

    function countdownSequence()
        cla(ax);
        % Get current axes limits
        xLim = get(ax, 'XLim');
        yLim = get(ax, 'YLim');
        centerX = mean(xLim);
        centerY = mean(yLim);

        countdown_text = {'Ready', 'Set', 'Go!'};
        for i = 1:length(countdown_text)
            cla(ax);
            text(centerX, centerY, countdown_text{i}, ...
                'Color', 'k', 'FontSize', 30, 'HorizontalAlignment', 'center', 'Parent', ax);
            drawnow;
            pause(1);
        end
    end

    function [score, isGameOver] = runGame()
        %% Clear previous graphics
        cla(ax);

        % Update axes limits before starting the game
        updateAxesLimits();

        % Get current axes limits
        xLim = get(ax, 'XLim');
        yLim = get(ax, 'YLim');

        %% Game Parameters
        electrode.pos = [xLim(1) + 20, 50];  % Starting position [x, y]
        electrode.vel = 0;           % Initial vertical velocity
        gravity = -0.1;              % Gravity effect (easier start)
        thrust = 0.4;                % Upward force when space bar is pressed
        max_vel = 5;                 % Maximum upward velocity
        obstacles = [];              % List to hold obstacles
        obstacle_speed = -0.5;       % Speed of obstacles moving left
        obstacle_gap = 40;           % Vertical gap in obstacles
        obstacle_freq = 200;         % Frequency of obstacle generation (frames)
        frame_rate = 30;             % Frames per second
        score = 0;                   % Player's score
        isGameOver = false;          % Game over flag
        frame_count = 0;             % Frame counter
        obstdist_count = 0;           % Obstacle distance counter
        difficulty_increase_interval = 1000; % Frames after which difficulty increases
        difficulty_level = 1;        % Initial difficulty level
        brain_shift_active = 0;      % brain shift active
        brain_shift_duration = 20;   % brain shift duration (frames)
        brain_shift_inited = 0;      % brain shift inited


        %% Set Key Functions for Game
        set(fig, 'KeyPressFcn', @keyDown, 'KeyReleaseFcn', @keyUp);

        %% Main Game Loop
        thisobst_dist=(obstacle_freq*rand);
        while ~isGameOver
            % Update frame count
            frame_count = frame_count + 1;
            obstdist_count = obstdist_count + 1;


            % Increase difficulty at intervals
            if mod(frame_count, difficulty_increase_interval) == 0
                difficulty_level = difficulty_level + 1;
                gravity = gravity - 0.03; % Increase gravity
                obstacle_speed = obstacle_speed - 0.3; % Increase obstacle speed
                if obstacle_gap > 20
                    obstacle_gap = obstacle_gap - 5; % Decrease gap size gradually
                end
                if obstacle_freq > 30
                    obstacle_freq = obstacle_freq - 30; % Increase obstacle frequency
                end

                % Notify user of level up
                notifyLevelUp(difficulty_level);
            end


            % Brain shift effect starts from level 4
            if difficulty_level >= 4
                if rand() < 0.001 || brain_shift_active % 5% chance of starting brain shift each frame
                    % bs init
                    brain_shift_active = 1;
                    if ~brain_shift_inited
                        brain_shift_start = frame_count;
                        brain_shift_inited = 1;
                    end
                    flashRectangle = rectangle('Position', [xLim(1), yLim(1), diff(xLim), diff(yLim)], ...
                        'FaceColor', [0.8, 0.8, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Parent', ax);

                    % Bring the flash rectangle to the front
                    uistack(flashRectangle, 'top');

                    xLim = get(ax, 'XLim');
                    yLim = get(ax, 'YLim');
                    centerX = mean(xLim);
                    bstext=text(centerX, yLim(2) - 20, 'Uh oh, brain shift!!', ...
                        'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Parent', ax);

                    % Update the figure to show the effect
                    drawnow;


                    % Shake electrode and obstacles during brain shift
                    shakeAmount = 2*randn(1, length(obstacles)); % Random shake magnitude
                    for i = 1:length(obstacles)
                        obstacles(i).pos = obstacles(i).pos + shakeAmount(i);
                    end

                    bstext.Position = bstext.Position + randn(1,3);

                end
            end

            % stop brain shift
            if brain_shift_active
                if frame_count-brain_shift_start >= brain_shift_duration
                    brain_shift_active = 0;
                    brain_shift_inited = 0;

                    % Remove the flash rectangle
                    delete(flashRectangle);
                    delete(bstext);

                end
            end

            % Handle zapping
            if zap_pressed && zaps_remaining > 0
                % Perform zapping action
                zaps_remaining = zaps_remaining - 1;
                % Call the function to handle zapping
                obstacles = zapObstacles(electrode, obstacles);
                % Reset zap_pressed to prevent continuous firing
                zap_pressed = false;
            elseif zap_pressed
                if ~isempty(splatSound)
                    sound(splatSound.y(1:5000), splatSound.Fs);
                end
            end

            % Handle input
            if space_pressed
                electrode.vel = min(electrode.vel + thrust, max_vel);
            else
                electrode.vel = electrode.vel + gravity;
            end

            % Update electrode position
            electrode.pos(2) = electrode.pos(2) + electrode.vel;

            % Keep electrode within vertical bounds
            if electrode.pos(2) < 0
                electrode.pos(2) = 0;
                electrode.vel = 0;
            elseif electrode.pos(2) > 100
                electrode.pos(2) = 100;
                electrode.vel = 0;
            end

            % Generate new obstacles
            if obstdist_count>thisobst_dist %mod(frame_count, obstacle_freq) == 0
                obstdist_count = 0; % reset
                thisobst_dist=(obstacle_freq*rand);
                gap_y = randi([30, 70]);
                neuron_size = obstacle_gap;

                obstacle_start_x = xLim(2) + 2;  % Start obstacles slightly off-screen

                % Top neuron obstacle
                top_neuron = struct('pos', [obstacle_start_x, gap_y + neuron_size / 2], ...
                    'size', [5, yLim(2) - (gap_y + neuron_size / 2)], ...
                    'type', 'neuron', 'subtype', 'top');

                % Bottom neuron obstacle
                bottom_neuron = struct('pos', [obstacle_start_x, 0], ...
                    'size', [5, gap_y - neuron_size / 2], ...
                    'type', 'neuron', 'subtype', 'bottom');

                obstacles = [obstacles, top_neuron, bottom_neuron];

                % Starting from level 2, add new obstacle types
                if difficulty_level >= 3
                    % Randomly decide whether to add blood vessel or glia
                    rand_choice = rand;
                    if rand_choice < 0.5
                        % Add a blood vessel
                        vessel_length = 100;  % Full width
                        vessel_thickness = 2 + rand * 2;  % Thickness between 2 and 4
                        vessel_y = randi([20, 80]);
                        blood_vessel = struct('pos', [100, vessel_y], ...
                            'size', [vessel_length, vessel_thickness], ...
                            'type', 'blood_vessel', 'subtype', 'none');
                        obstacles = [obstacles, blood_vessel];
                    else
                        % Add a glia cell
                        glia_radius = 5 + rand * 5;  % Radius between 5 and 10
                        glia_x = 100;
                        glia_y = randi([20, 80]);
                        glia_cell = struct('pos', [glia_x, glia_y], ...
                            'size', [glia_radius, glia_radius], ...
                            'type', 'glia', 'subtype', 'none');
                        obstacles = [obstacles, glia_cell];
                    end
                elseif difficulty_level == 2
                    % Add a glia cell
                    glia_radius = 5 + rand * 5;  % Radius between 5 and 10
                    glia_x = 100;
                    glia_y = randi([20, 80]);
                    glia_cell = struct('pos', [glia_x, glia_y], ...
                        'size', [glia_radius, glia_radius], ...
                        'type', 'glia', 'subtype', 'none');
                    obstacles = [obstacles, glia_cell];
                end
            end

            % Update obstacle positions
            for i = 1:length(obstacles)
                obstacles(i).pos(1) = obstacles(i).pos(1) + obstacle_speed;
            end

            % Remove off-screen obstacles
            obstacles = obstacles(arrayfun(@(obs) obs.pos(1) > (xLim(1) - 10), obstacles));

            % Collision detection
            for i = 1:length(obstacles)
                obs = obstacles(i);
                if checkCollision(electrode, obs)
                    isGameOver = true;
                    break;
                end
            end

            % Clear axes
            cla(ax);

            % Draw electrode (DBS electrode representation)
            drawElectrode(electrode.pos);

            % Draw obstacles
            for i = 1:length(obstacles)
                obs = obstacles(i);
                switch obs.type
                    case 'neuron'
                        drawNeuron(obs.pos, obs.size, obs.subtype);  % Now matches the function definition
                    case 'blood_vessel'
                        drawBloodVessel(obs);
                    case 'glia'
                        drawGlia(obs);
                end
            end

            % Update score
            score = score + 1;
            
            % Replenish zaps every 1000 points
            if mod(score, 600) == 0
                zaps_remaining = zaps_remaining + 1;
                % Display a notification
                text(50, 90, '+1 Zap!', 'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Parent', ax);
                drawnow;
                sound(gongSound.y(1:5000), gongSound.Fs);
                pause(0.5);  % Pause to let the player see the message
            end

            h = title(ax, ['Fiber Score: ', num2str(score), '   Level: ', num2str(difficulty_level), '   Zaps Remaining (.): ', num2str(zaps_remaining)], ...
                'Color', 'k', 'FontSize', 14);           
            % Pause for frame rate
            pause(1/frame_rate);

            % Check if figure is closed
            if ~ishandle(fig)
                isGameOver = true;
                break;
            end

            % Check if game is over (restart or exit)
            if isGameOver
                delete(h);
                break;
            end
        end
    end


    function obstacles = zapObstacles(electrode, obstacles)
        % Remove obstacles within a certain radius from the electrode
        zap_range = 40;  % Define the effective range of the zap

        % Loop through obstacles to check distances
        obstacles_to_remove = false(1, length(obstacles));  % Logical array to mark obstacles for removal

        for i = 1:length(obstacles)
            obs = obstacles(i);
            % Calculate distance between electrode and obstacle
            obs_pos = obs.pos;
            distance = sqrt((electrode.pos(1) - obs_pos(1))^2 + (electrode.pos(2) - obs_pos(2))^2);

            if distance <= zap_range
                % Mark obstacle for removal
                obstacles_to_remove(i) = true;

                % Optionally, play a sound or visual effect here
                % e.g., plot a flash or play a sound
            end
        end


        % Draw a semi-transparent rectangle over the axes to simulate the flash
        xLim = get(ax, 'XLim');
        yLim = get(ax, 'YLim');
        flashRectangle = rectangle('Position', [xLim(1), yLim(1), diff(xLim), diff(yLim)], ...
            'FaceColor', [1, 1, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Parent', ax);

        % Bring the flash rectangle to the front
        uistack(flashRectangle, 'top');

        % Draw electric effect
        drawElectricEffect(ax, electrode.pos);

        % Update the figure to show the effect
        drawnow;

        % Pause briefly to let the player see the effect
        pause(0.05);

        % Remove the flash rectangle
        delete(flashRectangle);

        if ~isempty(handelSound)
            sound(handelSound.y(1:5000), handelSound.Fs);
        end

        % Remove the marked obstacles
        obstacles(obstacles_to_remove) = [];

        % Optionally, you can add a visual effect to indicate zapping
    end

    function drawElectricEffect(ax, electrode_pos)
        % Draw an electric effect emanating from the electrode
        numBolts = 5;  % Number of lightning bolts
        maxBoltLength = 55;
        for i = 1:numBolts
            % Random angle for each bolt
            angle = rand * 2 * pi;
            % Random length between 50% and 100% of maxBoltLength
            length = maxBoltLength * (0.5 + 0.5 * rand);
            % Number of segments in the bolt
            numSegments = 5;
            % Initialize bolt coordinates
            x = electrode_pos(1);
            y = electrode_pos(2);
            boltX = x;
            boltY = y;
            for j = 1:numSegments
                % Random small angle deviation for zig-zag effect
                angleDeviation = (rand - 0.5) * pi / 6;
                segmentLength = length / numSegments;
                angle = angle + angleDeviation;
                x = x + segmentLength * cos(angle);
                y = y + segmentLength * sin(angle);
                boltX = [boltX, x];
                boltY = [boltY, y];
            end
            % Draw the bolt
            plot(ax, boltX, boltY, 'Color', [0, 0.5, 1], 'LineWidth', 2);
        end
    end

    function keyDown(~, event)
        if strcmp(event.Key, 'space')
            space_pressed = true;
        elseif strcmp(event.Key, 'period') || strcmp(event.Character, '.')
            zap_pressed = true;
        elseif strcmp(event.Key, 'escape')
            % Exit the game
            delete(fig);
        end
    end

    function keyUp(~, event)
        if strcmp(event.Key, 'space')
            space_pressed = false;
        elseif strcmp(event.Key, 'period') || strcmp(event.Character, '.')
            zap_pressed = false;
        end
    end

    function startKeyDown(~, event)
        if strcmp(event.Key, 'space')
            waitForStart = false;
        elseif strcmp(event.Key, 'escape')
            delete(fig);
        end
    end

    function restartKeyDown(~, event)
        if strcmp(event.Key, 'space')
            waitForRestart = false;
            isGameOver = false;
            space_pressed = false; % Reset space_pressed
            cla(ax);
        elseif strcmp(event.Key, 'escape')
            delete(fig);
        end
    end

    function closeGame(~, ~)
        stop(player);  % Stop the background music
        delete(fig);
    end

    function collision = checkCollision(electrode, obs)
        % Check if the electrode collides with an obstacle

        % Electrode dimensions
        electrode_width = 1.5;
        electrode_height = 6;

        % Electrode bounding box
        electrode_left = electrode.pos(1) - electrode_width / 2;
        electrode_right = electrode.pos(1) + electrode_width / 2;
        electrode_bottom = electrode.pos(2) - electrode_height / 2;
        electrode_top = electrode.pos(2) + electrode_height / 2;

        collision = false;

        switch obs.type
            case 'neuron'
                % Existing neuron collision logic
                x = obs.pos(1);
                y = obs.pos(2);
                scale = obs.size(2) / 50;
                soma_radius = 3 * scale;

                % Check collision with soma
                collision_with_soma = rectCircleCollision(electrode_left, electrode_right, electrode_bottom, electrode_top, ...
                    x, y, soma_radius);

                % Define axon bounding box
                if strcmp(obs.subtype, 'top')
                    axon_top = 100;
                    axon_bottom = y + soma_radius;
                else
                    axon_top = y - soma_radius;
                    axon_bottom = 0;
                end
                axon_left = x - 1;
                axon_right = x + 1;

                % Check collision with axon
                collision_with_axon = rectRectCollision(electrode_left, electrode_right, electrode_bottom, electrode_top, ...
                    axon_left, axon_right, axon_bottom, axon_top);

                collision = collision_with_soma || collision_with_axon;

            case 'blood_vessel'
                % Blood vessel collision detection
                vessel_left = obs.pos(1);
                vessel_right = obs.pos(1) + obs.size(1);
                vessel_thickness = obs.size(2);
                vessel_y = obs.pos(2);
                vessel_top = vessel_y + vessel_thickness / 2;
                vessel_bottom = vessel_y - vessel_thickness / 2;

                collision = rectRectCollision(electrode_left, electrode_right, electrode_bottom, electrode_top, ...
                    vessel_left, vessel_right, vessel_bottom, vessel_top);

            case 'glia'
                % Glia collision detection
                x = obs.pos(1);
                y = obs.pos(2);
                radius = obs.size(1);

                collision = rectCircleCollision(electrode_left, electrode_right, electrode_bottom, electrode_top, ...
                    x, y, radius);
        end
    end

    function collision = rectCircleCollision(rect_left, rect_right, rect_bottom, rect_top, circle_x, circle_y, circle_radius)
        % Check if a rectangle and a circle collide
        % Rectangle defined by left, right, bottom, top
        % Circle defined by center (circle_x, circle_y) and radius

        % Find the closest point on the rectangle to the circle center
        closest_x = max(rect_left, min(circle_x, rect_right));
        closest_y = max(rect_bottom, min(circle_y, rect_top));

        % Calculate the distance between the circle's center and this closest point
        distance_x = circle_x - closest_x;
        distance_y = circle_y - closest_y;

        % If the distance is less than the circle's radius, there's a collision
        distance_squared = (distance_x)^2 + (distance_y)^2;
        collision = distance_squared < (circle_radius)^2;
    end
    function collision = rectRectCollision(rect1_left, rect1_right, rect1_bottom, rect1_top, ...
            rect2_left, rect2_right, rect2_bottom, rect2_top)
        % Check if two rectangles collide
        collision = ~(rect1_right < rect2_left || ...
            rect1_left > rect2_right || ...
            rect1_top < rect2_bottom || ...
            rect1_bottom > rect2_top);
    end

    function drawElectrode(position)
        % Draw a DBS electrode at the given position
        x = position(1);
        y = position(2);

        % Electrode shaft
        %plot(ax, [x, x], [y - 5, y + 5], 'k-', 'LineWidth', 2);

        % Electrode contacts (4 rounded squares with equal sides)
        contact_positions = linspace(y - 3, y + 3, 4);
        side_length = 1.5;
        for cp = contact_positions
            rectangle('Position', [x - side_length/2, cp - side_length/2, side_length, side_length], ...
                'Curvature', 0.2, 'FaceColor', [0.2, 0.2, 0.2], 'EdgeColor', 'k', 'Parent', ax);
        end
    end

    function msc=getMusic(player)
        if player.CurrentSample<5629671
            msc='Dandy Junior / DNDY';
        else
            msc='Dandy Junior / DNDY & Neuston';
        end
    end


    function showMusicCredits
        try
            mscCredits=text(centerX, 40, ['Music by ',getMusic(player),'.'], ...
                'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Parent', ax);
            drawnow;
        end
    end

    function delMusicCredits
        try
            delete(mscCredits);
        end
    end

    function drawNeuron(pos, size, type)
        % Draw a neuron-like obstacle with a soma, dendrites, and a vertical axon extending to the edge

        x = pos(1);
        y = pos(2);
        scale = size(2) / 50;  % Scale the neuron size based on obstacle size

        % Parameters for neuron appearance
        soma_radius = 3 * scale;        % Radius of the cell body
        num_dendrites = 5;              % Number of dendrites
        dendrite_length = 10 * scale;   % Length of dendrites

        % Draw the soma (cell body)
        theta = linspace(0, 2*pi, 100);
        soma_x = x + soma_radius * cos(theta);
        soma_y = y + soma_radius * sin(theta);
        fill(ax, soma_x, soma_y, [0.8, 0.8, 0.8], 'EdgeColor', 'k');

        % Draw dendrites
        for i = 1:num_dendrites
            angle = rand * 2 * pi;  % Random angle for each dendrite
            dx = dendrite_length * cos(angle);
            dy = dendrite_length * sin(angle);

            % Dendrite branches
            num_branches = randi([1, 3]);
            branch_angles = angle + linspace(-pi/4, pi/4, num_branches);

            for ba = branch_angles
                bdx = (dendrite_length * 0.7) * cos(ba);
                bdy = (dendrite_length * 0.7) * sin(ba);

                % Draw the branch
                plot(ax, [x + dx, x + dx + bdx], [y + dy, y + dy + bdy], 'k-', 'LineWidth', 1);
            end

            % Draw the dendrite
            plot(ax, [x, x + dx], [y, y + dy], 'k-', 'LineWidth', 2);
        end

        % Draw the axon extending to the top or bottom
        if strcmp(type, 'top')
            % Neuron is at the top, axon extends upward
            axon_x = x;
            axon_y_start = y + soma_radius;
            axon_y_end = 100;  % Top of the screen
        else
            % Neuron is at the bottom, axon extends downward
            axon_x = x;
            axon_y_start = y - soma_radius;
            axon_y_end = 0;    % Bottom of the screen
        end
        plot(ax, [axon_x, axon_x], [axon_y_start, axon_y_end], 'k-', 'LineWidth', 2);
    end

    function drawBloodVessel(obs)
        % Draw a blood vessel obstacle

        x_start = obs.pos(1);
        y_start = obs.pos(2);
        length = obs.size(1);
        thickness = obs.size(2);

        % Generate a wavy line to represent the vessel
        num_points = 50;
        x_values = linspace(x_start, x_start + length, num_points);
        y_values = y_start + thickness * sin((x_values - x_start) / length * 4 * pi);  % 4 full waves

        % Draw the vessel
        plot(ax, x_values, y_values, 'r-', 'LineWidth', thickness);
    end

    function drawGlia(obs)
        % Draw a glia cell obstacle

        x = obs.pos(1);
        y = obs.pos(2);
        radius = obs.size(1);

        % Draw the main body
        theta = linspace(0, 2*pi, 100);
        glia_x = x + radius * cos(theta);
        glia_y = y + radius * sin(theta);
        fill(ax, glia_x, glia_y, [0.6, 0.8, 1], 'EdgeColor', 'k');  % Light blue color

        % Draw protrusions (processes)
        num_processes = 8;
        for i = 1:num_processes
            angle = (i / num_processes) * 2 * pi + (rand - 0.5) * 0.4;  % Slight randomness
            process_length = radius * (1 + rand * 0.5);
            dx = process_length * cos(angle);
            dy = process_length * sin(angle);
            plot(ax, [x + radius * cos(angle), x + dx], [y + radius * sin(angle), y + dy], 'k-', 'LineWidth', 1);
        end
    end

    function notifyLevelUp(level)
        % Function to notify the user of a level up

        % Play a pleasant sound
        chirp=load(fullfile(matlabroot, 'toolbox/matlab/audiovideo/chirp.mat'));
        S = warning('off', 'MATLAB:audiovideo:audioplayer:noAudioOutputDevice');
        sound(chirp.y(1000:-1:1),chirp.Fs/2);
        warning(S);

        xLim = get(ax, 'XLim');
        yLim = get(ax, 'YLim');
        centerX = mean(xLim);
        text(centerX, yLim(2) - 20, ['Level ', num2str(level)], ...
            'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Parent', ax);
        drawnow;
        pause(1);
    end
    function displayEndCredits(commitMessages)
        % Display the scrolling end credits with commit data

        % Combine commit messages into a single string
        creditText = strjoin(commitMessages, '\n\n');

        % Check if the text length is still too long
        maxTextLength = 10000;  % Maximum allowed text length
        if length(creditText) > maxTextLength
            creditText = creditText(1:maxTextLength);
            creditText = [creditText, '\n\n...'];
        end

        % Initial position of the text
        yPos = -50;  % Start below the visible area
        xPos = xLim(2)-(xLim(2)-xLim(1))/10; % 90%

        xLim = get(ax, 'XLim');

        % Create the text object
        creditHandle = text(xPos, yPos, creditText, 'Color', 'k', 'FontSize', 10, ...
            'HorizontalAlignment', 'right', 'Parent', ax, ...
            'Interpreter', 'none');

        % Get the extent (height) of the text
        drawnow;  % Update graphics to get accurate extent
        textExtent = get(creditHandle, 'Extent');
        textHeight = textExtent(4);

        % Scroll the text upwards
        % Calculate scroll parameters
        totalDistance = textHeight + 100;  % Height of text plus screen height
        desiredSpeed = 20;                 % Scroll speed in units per second (adjust as needed)
        totalTime = totalDistance / desiredSpeed;  % Total scrolling time in seconds
        framesPerSecond = 30;              % Target frame rate
        totalFrames = totalTime * framesPerSecond;
        scrollSpeed = totalDistance / totalFrames;  % Units per frame

        % Start scrolling
        while yPos < textHeight + 10 && waitForRestart  % Scroll until the text moves above the visible area or interrupted
            showMusicCredits;
            % Check if figure or text object is closed
            if ~ishandle(fig) || ~ishandle(creditHandle)
                return;
            end

            yPos = yPos + scrollSpeed;
            set(creditHandle, 'Position', [xPos, yPos, 0]);
            drawnow;

            pause(1 / framesPerSecond);  % Control frame rate
            delMusicCredits;
        end

    end

    function commitMessages = fetchCommitData()
        % Fetch the latest commits from the Lead-DBS GitHub repository

        % GitHub API base URL for commits
        baseUrl = 'https://api.github.com/repos/netstim/leaddbs/commits';

        % Parameters
        branch = 'develop';
        per_page = 100;  % Maximum per GitHub API
        max_commits = 100;  % Limit to 100 commits to prevent text overflow
        total_pages = ceil(max_commits / per_page);

        commitMessages = {};
        options = weboptions('Timeout', 15);

        % Loop over pages to fetch commits
        for page = 1:total_pages
            % Construct the API URL with pagination
            apiUrl = sprintf('%s?sha=%s&per_page=%d&page=%d', baseUrl, branch, per_page, page);

            try
                % Fetch data from GitHub API
                data = webread(apiUrl, options);

                % Extract commit messages
                numCommits = length(data);
                for i = 1:numCommits
                    commit = data(i).commit;
                    author = commit.author.name;
                    date = commit.author.date(1:10);  % Extract date only
                    message = commit.message;

                    % Truncate message to 200 characters
                    maxMessageLength = 200;
                    if length(message) > maxMessageLength
                        message = [message(1:maxMessageLength), '...'];
                    end

                    % Combine author, date, and message
                    commitMessages{end+1} = sprintf('Author: %s\nDate: %s\nMessage: %s\n', author, date, message);
                end

                % Check if we've reached the maximum number of commits
                if length(commitMessages) >= max_commits
                    commitMessages = commitMessages(1:max_commits);
                    break;
                end

            catch ME
                % Handle errors (e.g., network issues)
                commitMessages = {'Unable to fetch commit data.'};
                disp(['Error fetching commit data: ', ME.message]);
                break;
            end

            % Pause to avoid hitting rate limits
            pause(1);  % Adjust as needed based on rate limits
        end
    end

    function updateAxesLimits()
        % Get axes position in pixels
        set(ax, 'Units', 'pixels');
        axPos = get(ax, 'Position');
        axWidth = axPos(3);
        axHeight = axPos(4);

        % Calculate aspect ratio based on axes size
        aspectRatio = axWidth / axHeight;

        % Set y-axis limits
        yLim = [0, 100];

        % Adjust x-axis limits based on aspect ratio
        xRange = (yLim(2) - yLim(1)) * aspectRatio;
        xLim = [0, xRange];

        % Update axes limits
        xlim(ax, xLim);
        ylim(ax, yLim);

        % Update DataAspectRatio to maintain proportions
        set(ax, 'DataAspectRatio', [1 1 1]);

        % Reset axes Units to normalized
        set(ax, 'Units', 'normalized');
    end

end
