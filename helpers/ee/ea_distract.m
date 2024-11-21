function ea_distract()
    % Lead-DBS Neural Navigator: A minimal Easter Egg game within Lead-DBS controlled by the space bar.
    % The player navigates a DBS electrode through neural pathways, avoiding obstacles.

    %% Global Game Variables
    space_pressed = false;   % Space bar press flag
    isGameOver = false;      % Game over flag
    waitForRestart = false;  % Restart flag
    fig = [];                % Figure handle

    %% Set Up Figure Window
    fig = figure('Name', 'Lead-DBS Neural Navigator', 'NumberTitle', 'off', ...
                 'MenuBar', 'none', 'ToolBar', 'none', 'Color', 'w', ...
                 'KeyPressFcn', @keyDown, 'KeyReleaseFcn', @keyUp, ...
                 'CloseRequestFcn', @closeGame, ...
                 'WindowState', 'maximized');  % Maximize the window on launch
    ax = axes('Parent', fig);
    axis(ax, [0 100 0 100]);
    axis(ax, 'manual');
    hold(ax, 'on');
    set(ax, 'Color', 'w', 'XColor', 'none', 'YColor', 'none', 'DataAspectRatio', [1 1 1]);



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

            cla(ax);
            text(50, 60, 'Game Over!', ...
                'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Parent', ax);
            text(50, 50, ['Fiber Score: ', num2str(score)], ...
                'Color', 'k', 'FontSize', 16, 'HorizontalAlignment', 'center', 'Parent', ax);
            text(50, 40, 'Press Space to Restart or Esc to Exit', ...
                'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Parent', ax);

            % Wait for space bar to restart or escape to exit
            waitForRestart = true;
            set(fig, 'KeyPressFcn', @restartKeyDown);
            while waitForRestart
                pause(0.1);
                if ~ishandle(fig)
                    return; % Exit if figure is closed
                end
            end
        else
            break; % Exit the game loop if figure is closed
        end
    end

    %% Nested Functions

    function introScreen()
        % Display the intro screen with the provided ASCII art
        cla(ax);
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
        text(50, 50, ascii_art, 'Color', 'k', 'FontName', 'FixedWidth', ...
            'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Parent', ax, 'Interpreter', 'none');
        text(50, 10, 'Press Space to Start and to control the Electrode!', ...
            'Color', 'k', 'FontSize', 14, 'HorizontalAlignment', 'center', 'Parent', ax, 'Interpreter', 'none');
        drawnow;
    end

    function countdownSequence()
        % Display countdown: Ready, Set, Go
        cla(ax);
        countdown_text = {'Ready', 'Set', 'Go!'};
        for i = 1:length(countdown_text)
            cla(ax);
            text(50, 50, countdown_text{i}, ...
                'Color', 'k', 'FontSize', 30, 'HorizontalAlignment', 'center', 'Parent', ax);
            drawnow;
            pause(1);
        end
    end

    function [score, isGameOver] = runGame()
        %% Game Parameters
        electrode.pos = [20, 50];    % Starting position [x, y]
        electrode.vel = 0;           % Initial vertical velocity
        gravity = -0.1;              % Gravity effect (easier start)
        thrust = 0.4;                % Upward force when space bar is pressed
        max_vel = 5;                 % Maximum upward velocity
        obstacles = [];              % List to hold obstacles
        obstacle_speed = -0.5;       % Speed of obstacles moving left
        obstacle_gap = 60;           % Vertical gap in obstacles
        obstacle_freq = 120;         % Frequency of obstacle generation (frames)
        frame_rate = 30;             % Frames per second
        score = 0;                   % Player's score
        isGameOver = false;          % Game over flag
        frame_count = 0;             % Frame counter
        difficulty_increase_interval = 1000; % Frames after which difficulty increases
        difficulty_level = 1;        % Initial difficulty level

        %% Clear previous graphics
        cla(ax);

        %% Set Key Functions for Game
        set(fig, 'KeyPressFcn', @keyDown, 'KeyReleaseFcn', @keyUp);

        %% Main Game Loop
        while ~isGameOver
            % Update frame count
            frame_count = frame_count + 1;

            % Increase difficulty at intervals
            if mod(frame_count, difficulty_increase_interval) == 0
                difficulty_level = difficulty_level + 1;
                gravity = gravity - 0.03; % Increase gravity
                obstacle_speed = obstacle_speed - 0.2; % Increase obstacle speed
                if obstacle_gap > 20
                    obstacle_gap = obstacle_gap - 5; % Decrease gap size
                end
                if obstacle_freq > 60
                    obstacle_freq = obstacle_freq - 5; % Increase obstacle frequency
                end

                % Notify user of level up
                notifyLevelUp(difficulty_level);
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
            if mod(frame_count, obstacle_freq) == 0
                gap_y = randi([30, 70]);
                top_size = max(100 - (gap_y + obstacle_gap/2), 0);
                bottom_size = max(gap_y - obstacle_gap/2, 0);
                top_obstacle = struct('pos', [100, gap_y + obstacle_gap/2], 'size', [5, top_size]);
                bottom_obstacle = struct('pos', [100, 0], 'size', [5, bottom_size]);
                obstacles = [obstacles, top_obstacle, bottom_obstacle];
            end

            % Update obstacle positions
            for i = 1:length(obstacles)
                obstacles(i).pos(1) = obstacles(i).pos(1) + obstacle_speed;
            end

            % Remove off-screen obstacles
            obstacles = obstacles(arrayfun(@(obs) obs.pos(1) > -10, obstacles));

            % Collision detection
            for i = 1:length(obstacles)
                obs = obstacles(i);
                if checkCollision(electrode.pos, obs)
                    isGameOver = true;
                    break;
                end
            end

            % Clear axes
            cla(ax);

            % Draw electrode (DBS electrode representation)
            drawElectrode(electrode.pos);

            % Draw obstacles (Neuron-like structures)
            for i = 1:length(obstacles)
                obs = obstacles(i);
                drawNeuron(obs.pos, obs.size);
            end

            % Update score
            score = score + 1;
            title(ax, ['Score: ', num2str(score), '   Level: ', num2str(difficulty_level)], 'Color', 'k', 'FontSize', 14);

            % Pause for frame rate
            pause(1/frame_rate);

            % Check if figure is closed
            if ~ishandle(fig)
                isGameOver = true;
                break;
            end

            % Check if game is over (restart or exit)
            if isGameOver
                break;
            end
        end
    end

    function keyDown(~, event)
        if strcmp(event.Key, 'space')
            space_pressed = true;
        elseif strcmp(event.Key, 'escape')
            % Exit the game
            delete(fig);
        end
    end

    function keyUp(~, event)
        if strcmp(event.Key, 'space')
            space_pressed = false;
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

    function collision = checkCollision(pos, obs)
        % Check if the electrode collides with an obstacle
        collision = pos(1) + 1 >= obs.pos(1) && pos(1) - 1 <= obs.pos(1) + obs.size(1) && ...
                    pos(2) + 5 >= obs.pos(2) && pos(2) - 5 <= obs.pos(2) + obs.size(2);
    end

    function drawElectrode(position)
        % Draw a DBS electrode at the given position
        x = position(1);
        y = position(2);

        % Electrode shaft
        plot(ax, [x, x], [y - 5, y + 5], 'k-', 'LineWidth', 2);

        % Electrode contacts (4 rounded squares with equal sides)
        contact_positions = linspace(y - 3, y + 3, 4);
        side_length = 1.5;
        for cp = contact_positions
            rectangle('Position', [x - side_length/2, cp - side_length/2, side_length, side_length], ...
                      'Curvature', 0.2, 'FaceColor', [0.2, 0.2, 0.2], 'EdgeColor', 'k', 'Parent', ax);
        end
    end

    function drawNeuron(pos, size)
        % Draw a neuron-like obstacle
        x = pos(1);
        y = pos(2);
        height = size(2);

        % Ensure height is non-negative
        height = max(height, 0);

        % Avoid drawing if height is zero
        if height <= 0
            return;
        end

        % Draw a branching structure to represent a neuron
        % Main stem
        plot(ax, [x, x], [y, y + height], 'k-', 'LineWidth', 2);

        % Branches
        num_branches = 3;
        branch_positions = linspace(y, y + height, num_branches + 2);
        for bp = branch_positions(2:end-1)
            branch_length = 5;
            plot(ax, [x, x + branch_length], [bp, bp + 5], 'k-', 'LineWidth', 1);
            plot(ax, [x, x - branch_length], [bp, bp + 5], 'k-', 'LineWidth', 1);
        end
    end

    function notifyLevelUp(level)
        % Function to notify the user of a level up

        % Play a pleasant sound
        chirp=load(fullfile(matlabroot, 'toolbox/matlab/audiovideo/chirp.mat'));
        S = warning('off', 'MATLAB:audiovideo:audioplayer:noAudioOutputDevice');
        sound(chirp.y(1000:-1:1),chirp.Fs/2);
        warning(S);

        % Display level up message
        text(50, 80, ['Level ', num2str(level)], ...
            'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Parent', ax);
        drawnow;

        % Pause briefly to allow player to see the message
        pause(1);
    end
end
