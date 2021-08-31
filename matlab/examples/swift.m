function out = swift( )
%SWIFT sequence simulation
% (Please see the manuscript for a detailed description of what is simulated.)
%
% The preset parameter values reproduce figures 3,4 from the manuscript.
%
% =========================================================================
% It is recommended, to check out *nonzero* gradient rotation as well.
% Here, a few reasonable settings, which hopefully work on your machine,
% despite the proof-of-concept implementation of the simulator:
%
% n_TR = 20
% G_rot = 2
% n_profiles = 15
% epsilon = 1e-4 (Depending on available RAM, you may try smaller values.)
% final_dyn = 0

% (all other parameters can remain as preset)
%
% The simulation shows a systematic decrease (and a general instability) of
% the signal with increasing distance from the isocenter (x == 0).
% This is not surprising, since we have to deal with a rotating
% *unbalanced* gradient, which generates a variable local (and average
% intra-voxel) phase. The faster the rotation, the larger the effect (and
% problem).
% 
% Two remarks:

% - The apparent signal decay in Fig. 7 of Idiyatullin et al., JMR 2008:267-273. 
% seems to support this finding, although it may just as well have a
% different origin.
%
% - Eddy currents due to *non-cartesian* trajectories should have a similar
% effect on the bSSFP signal. While the residual gradient moment due to
% eddy currents is much smaller than in for SWIFT, the speed of rotation is
% usually considerably higher.
% A typical example would be MR fingerprinting with golden angle rotation
% of the spiral trajectories. (Whether the effect on the bSSFP signal is
% strong enough to be observable depends on the actual numbers though.)
% =========================================================================


% return value (currently not used)

out = [];

% use latex as default interpreter

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% field values

par = [];
opt = [];
str = [];

% tissue parameters

par.T1 = 20;
par.T2 = 10;
opt.T1 = [];
opt.T2 = [];
str.T1 = 'longitudinal relaxation time [ms]';
str.T2 = 'transverse relaxation time [ms]';

% overall sequence timing

par.TR = 1;
par.n_TR = 50; % 20; % recommended for rotating gradients;
par.ramp = 0.015;
opt.TR = [];
opt.n_TR = [];
opt.ramp = [];
str.TR = 'repetition time (== RF pulse duration) [ms]';
str.n_TR = 'number of RF pulses';
str.ramp = 'flip angle ramp to accelerate approach to steady-state (fraction of total simulation)';

% gradient rotation per TR

par.G_rot = 0; % 2;
opt.G_rot = [];
str.G_rot = 'Gradient rotation per TR [deg]';

% RF pulse parameters

par.fa = 10;
par.p_inc = 0;
par.pd_inc = 0;
opt.fa = [];
opt.p_inc = [];
opt.pd_inc = [];
str.fa = 'flip angle (per pulse) [deg]';
str.p_inc = 'phase increment (== neg. off-resonance) [deg]';
str.pd_inc = 'phase difference increment (RF spoiling) [deg]';

% hyperbolic secant (HS1) parameters

par.bw = 100;
opt.bw = [];
str.bw = 'HS1 pulse bandwidth [kHz] (full width)';

% displayed profiles and configurations in output

par.dx = 1000;
par.n_profiles = 2.5; % 15;
par.n_x = 501;
par.conf_min = -2;
par.conf_max = 2;
opt.dx = [];
opt.n_profiles = [];
opt.n_x = [];
opt.conf_min = [];
opt.conf_max = [];
str.dx = 'resolution of reconstructed profile [um]';
str.n_profiles = 'number of profiles';
str.n_x = 'number of spatial locations to reconstruct';
str.conf_min = 'smallest displayed bSSFP configuration';
str.conf_max = 'largest displayed bSSFP configuration';

% options for configuration model

par.alloc = 1e5;
par.epsilon = 1e-5; % 1e-4;
par.verbose = 'False';
opt.alloc = [];
opt.epsilon = [];
opt.verbose = { 'True', 'False' };
str.alloc = 'allocated space (will be increased, if necessary) [# of configurations]';
str.epsilon = 'discard configurations with L2 norm smaller than this (0 == max. accuracy)';
str.verbose = 'provide some informal output (stored and deleted configurations)';

% output every dwell time can be restricted

par.init_dyn = 0;
par.final_dyn = 1; % 0;
opt.init_dyn = [];
opt.final_dyn = [];
str.init_dyn = 'show intra-pulse dynamics at the beginning [# of TR]';
str.final_dyn = 'show intra-pulse dynamics at the end [# of TR]';

% conversion factor [deg] to [rad]

deg_2_rad = pi / 180;

% number of RF pulse segments

N_seg = 60;

% generate figure to show progress during the simulation

fig_bSSFP = figure;
fig_occ = figure;

% program loop

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    %% for simplicity, we assume the RF pulse duration to be equal to TR
    
    par.Tp = par.TR;
    
    %% check and set preparation ramp
    
    if ( par.ramp < 0 || par.ramp > 1 )
        
        fprintf( 1, 'ramp should be in interval [0,1], please correct.' );
        continue;
        
    end
    
    % flip angle ramp
    
    ramp = ones( 1, par.n_TR );
    
    if ( par.ramp > 0 )
        
        n_ramp = ceil( par.ramp * par.n_TR );
        
        ramp( 1 : n_ramp ) = sin( linspace( 1 / n_ramp, 1, n_ramp ) .* ( 0.5 * pi ) );
        
    end
    
    %% round TR to multiple of dwell time
    
    dw = par.Tp / N_seg;           % dwell time
    par.TR = N_seg * dw;               % round TR, if necessary

    %% check for bandwidth    
    % the Nyquist condition generates an upper bound for the whole FOV bandwidth
    
    bw_fov = 1 / dw;
    
    % We require the RF pulse frequency sweep to be at least that large
    
    if ( bw_fov > par.bw )
        
        fprintf( 1, '\n*****************************************************\n' );
        fprintf( 1, 'bw_fov = 1 / dwell time = %.2f larger than frequency sweep\n', bw_fov );
        fprintf( 1, 'Increase bw, Tp or reduce N_seg.' );
        fprintf( 1, '\n*****************************************************\n\n' );
        continue;
        
    else
        
        fprintf( 1, 'bw_fov < bw: %.1f < %.1f  --  OK\n', bw_fov, par.bw );
        
    end

    %% generate HSn pulse and check for rapid linear region
    
    HSn_par = par;
    HSn_par.n = 1;         % we only simulate HS1 pulses
    HSn_par.N_seg = N_seg; % number of RF pulse segments
    HSn_par.L_over = 1;    % no gapped pulse
    HSn_par.d_c = 1;       % full duty cycle
    
    HSn_pulse = HSn( HSn_par );
    
    % rapid passage condition: a / omega_1^2 << 1
    % a = sweep rate, calculated at center of RF pulse
    
    fprintf( 1, 'rapid linear region: %f >> 1 ?\n', HSn_pulse.a / HSn_pulse.omega_1_max^2 );
    
    % parameters needed for the output
    
    k_space = pi / par.dx;
    p_min = 2 * k_space * par.conf_min;
    p_max = 2 * k_space * par.conf_max;
    tau_min = dw * round( N_seg * par.conf_min );
    tau_max = dw * round( N_seg * par.conf_max );
    N_tau = round( ( tau_max - tau_min ) / dw ) + 1;
    p_rng = linspace( p_min, p_max, N_tau );
    tau_rng = linspace( tau_min, tau_max, N_tau );
    t_rng = linspace( 0, par.TR - dw, N_seg );
    
    % some quantities only need to be evaluated once
    
    t_exc_fov = [];             % effective excitation time with respect to FOV bandwidth
    bSSFP_profile = [];         % approximate steady-state bSSFP profile
    bSSFP_config = [];          % approximate steady-state bSSFP configuration
    
    %% initialize configuration model (SWIFT and pulsed FT NMR)
    
    cm_swift = CoMo;
    cm_dess = CoMo;
    
    % tissue parameters
    
    cm_swift.R1 = 1 / par.T1;
    cm_swift.R2 = 1 / par.T2;
    cm_swift.D = 0;
    
    cm_dess.R1 = 1 / par.T1;
    cm_dess.R2 = 1 / par.T2;
    cm_dess.D = 0;
    
    % allocated space
    
    cm_swift.alloc = round( par.alloc );
    cm_dess.alloc = round( par.alloc );
    
    % gradient moment per dwell time
    
    p_tot = 2 * pi / par.dx;  % net gradient moment per TR
    
    dp_dw_swift = p_tot / N_seg;
    
    dp_dw_dess = p_tot / ( N_seg / 3 );
    
    % for rotating gradients, we want a finer grid in configuration space
    
    if ( par.G_rot == 0 ) % 1d case
        
        dp_cm_swift = dp_dw_swift;             % same grid in CM
        dp_cm_dess = dp_dw_dess;             % same grid in CM
        
    else  % rotating gradients == 2d p-space
        
        % in presence of rotating gradients, the simulated resolution of
        % configuration space should be better than the increment per dwell
        % time

        dp_res = 100;
                               
        dp_cm_swift = dp_dw_swift / dp_res;
        dp_cm_dess = dp_dw_swift / dp_res;
        
    end
    
    % apply settings to CM
    
    cm_swift.d_tau = dw;
    cm_dess.d_tau = dw;
    
    if ( par.G_rot == 0 ) % 1d case
        
        cm_swift.d_p = [ dp_cm_swift; 0; 0 ];
        cm_dess.d_p = [ dp_cm_dess; 0; 0 ];
        
    else  % rotating gradients == 2d p-space
        
        cm_swift.d_p = 0.05 .* [ dp_cm_swift; dp_cm_swift; 0 ];
        cm_dess.d_p = 0.05 .* [ dp_cm_dess; dp_cm_swift; 0 ];
        
    end
        
    % accuracy, verbosity
    
    cm_swift.epsilon = par.epsilon;
    cm_dess.epsilon = par.epsilon;
    
    if ( isequal( par.verbose, 'True' ) )
        
        cm_swift.verbose = true;
        cm_dess.verbose = true;
        
    else
        
        cm_dess.verbose = false;
        
    end
    
    % invariant RF pulse parameters
    
    RF_par_swift = [];
    RF_par_dess = [];
    
    RF_par_dess.FlipAngle = par.fa * deg_2_rad;
    j_RF_dess = round( N_seg / 2 + 1 );   % where is the RF pulse located?
    dp_dess = zeros( N_seg, 1 );
    dp_dess( 1 : 20 ) = dp_dw_dess;  
    dp_dess( 21 : 40 ) = - dp_dw_dess;
    dp_dess( 41 : 60 ) = dp_dw_dess;
      
    % invariant time interval parameters
    
    % small time intervals, related to transmit periods
    
    Time_B1_swift = [];
    Time_B1_swift.tau = cm_swift.d_tau;
    
    Time_dess = [];
    Time_dess.tau = cm_dess.d_tau;
    
    % initialize phase increment
    
    phase_inc_rad = 0;       % actual phase increment
    
    p_inc_rad = par.p_inc * deg_2_rad;
    pd_inc_rad = par.pd_inc * deg_2_rad;
    
    % gradient rotation
    % different from SWIFT, we rotate the gradient contiuously, i.e. every dwell time
    % to allow for TR = Tp
    % (should not make much difference for TR > Tp and
    % generates an even more silent sequence..)

    phiG = ( ( - par.n_TR * N_seg + 1 ) : 0 ) .* ( par.G_rot / N_seg );
    cG = cosd( phiG );
    sG = sind( phiG );

    %% execute sequence
    % start timing
    
    tic;
    
    % allocate space for the configurations
    
    N_dwell = par.n_TR * N_seg; % total number of intervals to simulate
    
    occ_config_swift = zeros( 2 * N_dwell + 1, N_dwell );
    occ_config_dess = zeros( 2 * ( N_dwell / 3 ) + 1, N_dwell );
        
    % reduced ranges for more detailed display
    
    occ_rng_det_swift = N_dwell + 1 + ( par.conf_min * N_seg : par.conf_max * N_seg );
    occ_rng_det_dess = ( N_dwell / 3 ) + 1 + ( par.conf_min * ( N_seg / 3 ) : par.conf_max * ( N_seg / 3 ) );

    n_detail = 3;
    time_rng_det = N_dwell - n_detail * N_seg : N_dwell;
    
    for j_TR = 1 : par.n_TR
        
        fprintf( 1, 'RF pulse %d / %d\n', j_TR, par.n_TR );        

        % first, we start with the RF pulse

        for j_seg = 1 : N_seg
            
            % set RF parameters (for each segment)
            
            RF_par_swift.FlipAngle = ramp( j_TR ) * HSn_pulse.alpha( j_seg );
            RF_par_swift.Phase = HSn_pulse.phi( j_seg ) + phase_inc_rad;
            
            % RF (sub)pulse
            
            cm_swift.RF( RF_par_swift );
            
            % the RF pulse of the DESS sequence shall be applied in the
            % center of the TR interval
            
            if ( j_seg == j_RF_dess )
                
                RF_par_dess.Phase = phase_inc_rad;
                
                cm_dess.RF( RF_par_dess );
                
            end
            
            % dwell time
            % the gradient moment may rotate
            
            j_tot = ( j_TR - 1 ) * N_seg + j_seg;
            
            Time_B1_swift.p = dp_dw_swift .* [ cG( j_tot ); sG( j_tot ); 0 ];
            
            cm_swift.time( Time_B1_swift );
            
            Time_dess.p = dp_dess( j_seg ) .* [ cG( j_tot ); sG( j_tot ); 0 ];
            
            cm_dess.time( Time_dess );
            
            %% output
            
            % for the manuscript, we want to display the configurations
            % every time point
            
            p_occ_swift = cm_swift.p( 1, cm_swift.occ );
            i_occ_swift = round( p_occ_swift / dp_dw_swift ) + N_dwell + 1;
            mp_occ_swift = sqrt( 2 ) .* cm_swift.m( 1, cm_swift.occ );
            
            occ_config_swift( i_occ_swift, j_tot ) = mp_occ_swift;

            b_dess = cm_dess.occ;
            b_dess( b_dess & cm_dess.m( 1, : ) == 0 ) = false;
            p_occ_dess = cm_dess.p( 1, b_dess );
            i_occ_dess = round( p_occ_dess / dp_dw_dess + ( N_dwell / 3 ) + 1 );
            mp_occ_dess = sqrt( 2 ) .* cm_dess.m( 1, b_dess );
            
            occ_config_dess( i_occ_dess, j_tot ) = mp_occ_dess;

            if ( j_TR <= par.init_dyn || j_TR > par.n_TR - par.final_dyn && j_seg ~= N_seg )
                
                    update_config = false;
                    compare_to_bSSFP( j_seg * dw );
                    
            end
            
            if ( j_seg == N_seg )
                
                if ( j_TR < par.n_TR )
                    
                    update_config = true;
                    
                else
                    
                    update_config = false;
                    
                end
                    
                compare_to_bSSFP( j_seg * dw );
                
            end
            
        end
        
        % update phase cycling
        
        phase_inc_rad = phase_inc_rad + p_inc_rad + j_TR * pd_inc_rad;
        
    end
    
    % stop timing
    
    toc;
    
end

%% helper functions to display results

    function compare_to_bSSFP ( t_ )
        
        % plots actual spatial bSSFP profile (left) and occupation in
        % configuration space (right)
        
        % the time t_ is a multiple of the dwell time dw and bound to the
        % interval [0,TR]
        
        %% comparison of local profile
        % spatial extent and locations
        
        dx_disp = par.dx * par.n_profiles;
        x_loc = linspace( - 0.5 * dx_disp, 0.5 * dx_disp, par.n_x );
        
        % CCM simulation ...
        
        sum_par = [];
        sum_par.x = [ cG( j_tot ); sG( j_tot ); 0 ] .* x_loc; % the profile rotates with the gradient
        
        res = cm_swift.sum( sum_par );
        mp_sim = squeeze( res.xy );
        
        % ... compared with bSSFP approximation
        % (for simplicity, absence of phase cycling is assumed)
        
        if ( isempty( bSSFP_profile ) )  % we need to calculate this only once
            
            fprintf( 1, 'Initialize bSSFP profile ...' );
            
            gamG = 2 * pi / ( par.TR * par.dx );  % gradient for 2 \pi dephasing over voxel per TR
            omega_fov = gamG .* x_loc;
            
            t_exc_fov = calc_t_exc( omega_fov );
            calc_bSSFP_profile( x_loc, 0, - t_exc_fov );

            fprintf( 1, ' done\n' );
            
        end
                        
        %% comparison of occupied configurations
        
        % for the configurations we distinguish between
        % 1d (constant gradients) and 2d (rotating gradients)

        x_loc_mm = x_loc ./ 1000;            

        if ( par.G_rot == 0 ) % constant gradient
            
            p_occ_swift = cm_swift.p( 1, cm_swift.occ );
            mp_occ_swift = sqrt( 2 ) .* cm_swift.m( 1, cm_swift.occ );
            
            b_discard = p_occ_swift < p_min | p_occ_swift > p_max;
            p_occ_swift( b_discard ) = [];
            mp_occ_swift( b_discard ) = [];
            
            % calculate theoretical bSSFP configuration
            
            if ( update_config && length( p_occ_swift ) == N_tau )
                
                calc_bSSFP_config( mp_occ_swift );
                
            end

            p_ = p_occ_swift ./ ( 2 * k_space );
            
            % output
            
            % show the plots (profiles vs. configurations) ...
            
            figure( fig_bSSFP );

            ro_ = 2;
            co_ = 2;
            
            y_min = 0;
            y_max_profile = 1.05 * max( max( abs( mp_sim ) ), max( abs( bSSFP_profile ) ) );
            y_max_conf = 1.05 * max( abs( mp_occ_swift ) );
            
            t_idx = round( mod( t_, par.TR ) / dw ) + 1;
            
            i_extra = 0;
            
            if ( t_idx == round( N_seg / 2 ) ) % show the center of the RF pulse separately

                i_extra = 2;
                
            end
                
            subplot( ro_, co_, 1 + i_extra );
            plot( x_loc_mm, abs( mp_sim ), x_loc_mm, abs( bSSFP_profile ) );
            ylim( [ y_min, y_max_profile ] );
            
            if ( i_extra == 0 )
                
                title( '$\left|\,m_+\left(x,0\right)\,\right|$' );

            else
            
                title( '$\left|\,m_+\left(x,\textrm{TR}/2\right)\,\right|$' );
                xlabel( '$x/\Delta x$' );

            end
            
            xlim( [ min( x_loc_mm ), max( x_loc_mm ) ] );
            xticks( - floor( par.n_profiles / 2 ) : floor( par.n_profiles / 2 ) );
                
            subplot( ro_, co_, 2 + i_extra );

            if ( isempty( bSSFP_config ) )
            
                plot( p_, abs( mp_occ_swift ) );

            else
                
                plot( p_, abs( mp_occ_swift ), p_rng ./ ( 2 * k_space ), abs( bSSFP_config( t_idx, : ) ) );
                
                if ( i_extra == 0 )
                    
                    title( '$\left|\,\hat{m}_+\left(p,0\right)\,\right|$' );
                    legend( 'CCM', 'bSSFP', 'Location', 'northwest' );
                    
                else
                    
                    title( '$\left|\,\hat{m}_+\left(p,\textrm{TR}/2\right)\,\right|$' );
                    xlabel( '$p\,\Delta x/2\pi$' );
                    
                end                
                
                xticks( par.conf_min : par.conf_max );
                
            end

            ylim( [ y_min, y_max_conf ] );

            drawnow;
            
            % ...and the images (occupied configurations)
            
            figure( fig_occ );

            ro_ = 2;
            co_ = 2;
            
            coma = colormap;
            coma( 1, : ) = 0;
            colormap( coma );

            col_max_swift = max( max( abs( occ_config_swift( occ_rng_det_swift, time_rng_det ) ) ) );
            col_max_dess = max( max( abs( occ_config_dess( occ_rng_det_dess, time_rng_det ) ) ) );
            
            if ( par.n_TR <= 10 )
                
                step_TR = 1;
                
            elseif ( par.n_TR <= 20 )
                
                step_TR = 5;
            
            else
                
                step_TR = 10;
                
            end
                
            % SWIFT all times *********************************************
            
            subplot( ro_, co_, 1 );            

            x = [ 0, par.n_TR ];
            y = [ - par.n_TR, par.n_TR ];
            x_det = [ par.n_TR - n_detail, par.n_TR ];
            y_det = [ par.conf_min, par.conf_max ];
            
            if ( col_max_swift )
                
                imagesc( x, y, abs( occ_config_swift ), [ 0, col_max_swift ] );
                
            else
                
                imagesc( x, y, abs( occ_config_swift ) );
                
            end

            title( '$\left|\,\hat{m}_+\left(p,t\right)\,\right|$' );
            set( gca, 'YDir', 'normal' );
            xticks( 0 : step_TR : par.n_TR );
            ylabel( '$p\,\Delta x/2\pi$' );

            % DESS all times *********************************************

            subplot( ro_, co_, 3 );
            
            if ( col_max_dess )
                
                imagesc( x, y, abs( occ_config_dess ), [ 0, col_max_dess ] );
                
            else
                
                imagesc( x, y, abs( occ_config_dess ) );
                
            end
                        
            set( gca, 'YDir', 'normal' );
            xticks( 0 : step_TR : par.n_TR );
            xlabel( '$t / \mathrm{TR}$' );
            ylabel( '$p\,\Delta x/2\pi$' );
            
            % SWIFT detail *********************************************
            
            subplot( ro_, co_, 2 );
            
            if( col_max_swift )
                
                imagesc( x_det, y_det, abs( occ_config_swift( occ_rng_det_swift, time_rng_det ) ), [ 0, col_max_swift ] );
                
            else
                
                imagesc( x_det, y_det, abs( occ_config_swift( occ_rng_det_swift, time_rng_det ) ) );
                
            end
            
            set( gca, 'YDir', 'normal' );
            xticks( par.n_TR - n_detail : par.n_TR );
            title( '$\left|\,\hat{m}_+\left(p,t\right)\,\right|$' );

            % DESS detail *********************************************
            
            subplot( ro_, co_, 4 );
            
            if ( col_max_dess )
                
                imagesc( x_det, y_det, abs( occ_config_dess( occ_rng_det_dess, time_rng_det ) ), [ 0, col_max_dess ] );
                
            else
                
                imagesc( x_det, y_det, abs( occ_config_dess( occ_rng_det_dess, time_rng_det ) ) );
                
            end
            
            set( gca, 'YDir', 'normal' );
            xlabel( '$t / \mathrm{TR}$' );
            xticks( par.n_TR - n_detail : par.n_TR );
            
            set( gcf, 'color', 'k' );
            set( gcf, 'InvertHardcopy', 'off' );

            drawnow;
            
        else % rotating gradient

            p_occ_swift = cm_swift.p( 1 : 2, cm_swift.occ );
            mp_occ_swift = sqrt( 2 ) .* cm_swift.m( 1, cm_swift.occ );
            
            b_discard = min( p_occ_swift ) < p_min | max( p_occ_swift ) > p_max;
            p_occ_swift( :, b_discard ) = [];
            mp_occ_swift( b_discard ) = [];
            
            % for rotating gradients, we do not want to compare with
            % a theoretical approximation, but just visually show the
            % occupation in configuration space
            
            n_2d = round( ( p_max - p_min ) / dp_dw_swift ) + 1;
            mp_2d = zeros( n_2d );

            ip_occ = round( ( p_occ_swift - p_min ) ./ ( p_max - p_min ) .* ( n_2d - 1 ) ) + 1;
            
            fprintf( 1, 'min( ip_occ ) = %d, max( ip_occ ) = %d\n', min( ip_occ( : ) ), max( ip_occ( : ) ) );
            
            for i = 1 : length( mp_occ_swift )
            
                mp_2d( sub2ind( [ n_2d, n_2d ], ip_occ( 2, i ), ip_occ( 1, i ) ) ) = ...
                    mp_2d( sub2ind( [ n_2d, n_2d ], ip_occ( 2, i ), ip_occ( 1, i ) ) ) + mp_occ_swift( i );
            
            end

            ro_ = 1;
            co_ = 1;
            
            figure( fig_bSSFP );

            y_min = 0;
            y_max_profile = 1.05 * max( max( abs( mp_sim ) ), max( abs( bSSFP_profile ) ) );
            
            subplot( ro_, co_, 1 );
            plot( x_loc_mm, abs( mp_sim ), x_loc_mm, abs( bSSFP_profile ) );
            ylim( [ y_min, y_max_profile ] );
            
            title( '$\left|\,m_+\left(x,0\right)\,\right|$' );
                
            xlabel( '$x/\Delta x$' );
                
            xlim( [ min( x_loc_mm ), max( x_loc_mm ) ] );
            xticks( - floor( par.n_profiles / 2 ) : floor( par.n_profiles / 2 ) );            
            
            figure( fig_occ );
            
            subplot( ro_, co_, 1 );
            
            coma = colormap;
            coma( 1, : ) = 0;
            colormap( coma );
            
            x_det = [ par.conf_min, par.conf_max ];
            y_det = [ par.conf_min, par.conf_max ];
            
            imagesc( x_det, y_det, abs( mp_2d ) );
            
            title( '$\left|\,\hat{m}_+\left({\bf p},t\right)\,\right|$' );
            set( gca, 'YDir', 'normal' );
            xlabel( '$p_y\,\Delta y/2\pi$' );
            ylabel( '$p_x\,\Delta x/2\pi$' );
            xticks( par.conf_min : par.conf_max );
            yticks( par.conf_min : par.conf_max );
            
            drawnow;
            
        end
        
    end

    function calc_bSSFP_profile ( x, om, ttx )
        %BSSFP_PROFILE return balanced SSFP
        
        % check arguments
        
        if ( length( x ) ~= length( ttx ) )
            
            error( 'Error in bSSFP_profile: x and ttx must have same length.' );
            
        end
        
        % adjust size
        
        n_x = length( x );
        
        x = reshape( x, [ n_x, 1 ] );
        ttx = reshape( ttx, [ n_x, 1 ] );
        
        % guarantee correct range of ttx
        
        ttx = mod( ttx, par.TR );
        
        % calculate the configurations according to
        % Ganter, Magn Reson Med 2006, 56:687-691
        
        E1 = exp( - par.TR / par.T1 );
        E2 = exp( - par.TR / par.T2 );
        c = cosd( par.fa );
        s = sind( par.fa );
        a = 1 - E1 * c;
        b = c - E1;
        E22 = E2^2;
        Lambda = ( a - b * E22 - sqrt( ( a^2 - E22 * b^2 ) * ( 1 - E22 ) ) ) / ( a - b );
        A = 0.5 * ( a + b ) * E2 / ( a - 0.5 * Lambda * ( a - b ) );
        out.n_max = round( - 4 / log10( A ) );  % to get a reasonably accurate (10^(-4)) result
        out.n = - out.n_max : out.n_max;
        idx_neg = 1 : out.n_max;
        out.m_n = A.^abs( out.n );
        out.m_n( idx_neg ) = - ( Lambda / ( E2 * A ) ) .* out.m_n( idx_neg );
        out.m_n = ( - 1i * ( 1 - E1 ) / ( a + b * Lambda ) * s ) .* out.m_n;
        
        % evaluate the sum at locations x
        
        gamG = 2 * pi / ( par.TR * par.dx );  % gradient for 2 pi dephasing over voxel per TR
        
        bSSFP_profile = exp( - ttx ./ par.T2 ) .* ...
            sum( ...
            exp( - 1i .* ( om + gamG .* x ) .* ( ttx + out.n .* par.TR ) ) .* out.m_n ...
            , 2 );
        
    end

    function calc_bSSFP_config( mp_occ_ )
        
        % set parameters
        
        A = pi * par.bw;
        beta = asech( 0.01 );

        R2 = 1 / par.T2;
       
        j_max = floor( par.bw * par.TR / 2 * tanh( beta ) );

        if ( abs( j_max - par.TR * par.bw / 2 ) < eps )
            
            j_max = j_max - 1;
            
        end

        % spending separate dimensions for independent variables

        t_rng = reshape( t_rng, [ length( t_rng ), 1 ] );
        tau_rng = reshape( tau_rng, [ 1, length( tau_rng ) ] );
        E2 = exp( - R2 .* ( mod( tau_rng, par.TR ) ) );

        % the configuration order depends on tau
            
        n_tau = round( ( tau_rng - mod( tau_rng, par.TR ) ) ./ par.TR );
        
        % search for best match
        
        j_opt = 0;
        norm_opt = 0;
        dev_min = Inf;
        
        for j_test = 1 : j_max
        
            j = - j_test : j_test;
            
            % spending separate dimensions for independent variables
            
            j = reshape( j, [ 1, 1, length( j ) ] );
            
            % location of poles and associated effective excitation times
            
            zj = ( 1i * 2 * pi / par.TR ) .* j - R2;
            tj = ( atanh( - 1i .* zj ./ A ) ./ beta + 1 ) .* ( 0.5 * par.Tp );
            
            % summing up residues
            
            ei_fak = sum( exp( ( 1i * 2 * pi / par.TR ) .* j .* ( tau_rng + tj ) ), 3 );
            
            % putting everything together
            
            res = ei_fak .* E2 .* out.m_n( n_tau + out.n_max + 1 ) ./ ( 2 * pi * par.TR );
            
            norm_fac = ( abs( res ) * abs( mp_occ_( : ) ) ) ./ sum( abs( res( : ) ).^ 2 );
            
            dev_ = sum( ( norm_fac .* abs( res( : ) ) - abs( mp_occ_( : ) ) ).^ 2 );
            
            if ( dev_ < dev_min )
                
                j_opt = j_test;
                dev_min = dev_;
                norm_opt = norm_fac;
                
            end

        end
                    
        fprintf( 'j_opt = %d, j_max = %d\n', j_opt, j_max );
        
        j = - j_opt : j_opt;
        
        j = reshape( j, [ 1, 1, length( j ) ] );
        
        % location of poles and associated effective excitation times
        
        zj = ( 1i * 2 * pi / par.TR ) .* j - R2;
        tj = ( atanh( - 1i .* zj ./ A ) ./ beta + 1 ) .* ( 0.5 * par.Tp );
        
        % summing up residues
        
        ei_fak = sum( exp( ( 1i * 2 * pi / par.TR ) .* j .* ( tau_rng - t_rng + tj ) ), 3 );
        
        % putting everything together
        
        bSSFP_config = norm_opt .* ei_fak .* E2 .* out.m_n( n_tau + out.n_max + 1 ) ./ ( 2 * pi * par.TR );
        
    end

    function t_exc_ = calc_t_exc ( omega_fov )
        % calculate effective excitation time
        
        %% effective excitation time
        % for these locations and the previous specification of the HSn
        % pulse, we can now estimate the position dependent excitation time,
        % based upon the local resonance condition
        % (which again depends on the gradient polarity)
        
        n_ = numel( omega_fov );
        
        t_exc_ = NaN( size( omega_fov ) );      % no excitation (NaN) outside the RF pulse bandwidth (should not happen)
        
        for j = 1 : n_
            
            i_lo = find( omega_fov( j ) > HSn_pulse.omega_RF, 1, 'last' );
            i_hi = find( omega_fov( j ) < HSn_pulse.omega_RF, 1 );
            
            if ( ~isempty( i_lo ) && ~isempty( i_hi ) )
                
                om_lo = HSn_pulse.omega_RF( i_lo );
                om_hi = HSn_pulse.omega_RF( i_hi );
                
                if ( om_lo > om_hi )
                    
                    error( 'Opposite orientation of frequency sweep. (Should not happen.)' );
                    
                end
                
                x_ = ( omega_fov( j ) - om_lo ) / ( om_hi - om_lo );  % linear interpolation
                
                t_exc_( j ) = ( 1 - x_ ) * HSn_pulse.t( i_lo ) + x_ * HSn_pulse.t( i_hi );
                
            end
            
        end
        
    end

end
    