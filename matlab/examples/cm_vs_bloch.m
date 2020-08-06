%% Comparison CoMoTk vs Bloch simulation for a random RF train
% 
% - no gradients, i.e. all magnetization pathways contribute
% - flip angle and phases randomized for each RF pulse
% - time intervals, randomly chosen from two different durations
% - single-peak tissue with freely adjustable off-resonance
% - no susceptibility effects (i.e. microscopic solution)
% - comparison with approximate solution (neglecting small configurations)
% - visualization of occupied vs stored configurations 

par = [];
opt = [];
str = [];

par.T1 = 100;
par.T2 = 10;
par.tau = [ 0.8, 1.2 ];
par.omega = 1;
par.n_inter = 500;
par.epsilon = 1e-6;

opt.T1 = [];
opt.T2 = [];
opt.tau = [];
opt.omega = [];
opt.n_inter = [];
opt.epsilon = [];

str.T1 = '[ms]';
str.T2 = '[ms]';
str.tau = '[ms] manual setting of *two* time interval durations';
str.omega = '[rad/ms] dephasing angle due to off-resonance (randomly chosen, if empty';
str.n_inter = 'number of RF pulses (separated by non-equidistant time intervals)';
str.epsilon = 'discard configuration vectors with L2 norm smaller than this';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
        
    if ( length( par.tau ) ~= 2 )
        
        fprintf( 1, 'length( tau ) == 2 required!' );
        
        continue;
           
    else
        
        tau = par.tau;
        
    end
        
    % off-resonance frequency and associated accumulated phase in the d intervals
    
    if ( isempty( par.omega ) )
    
        % random setting, if not specified
        
        omega = 2 * pi * rand;
        
    else
        
        % otherwise set by the user
        
        omega = par.omega;
        
    end
    
    %% initialize Bloch simulations
    
    % dimension of configuration model
    
    d = 2;
        
    % dephasing angle due to off-resonance, depending on the d durations tau(:)
    
    theta = - omega .* tau;
    
    % d relaxation and precession precession matrices
    
    E = zeros( 3, 3, d );
    
    E1 = exp( - tau / par.T1 );
    E2 = exp( - tau / par.T2 );
    
    E( 1, 1, : ) = E2;
    E( 2, 2, : ) = E2;
    E( 3, 3, : ) = E1;
    
    R_z = zeros( 3, 3, d );
    
    R_z( 1, 1, : ) = cos( theta );
    R_z( 2, 2, : ) = R_z( 1, 1, : );
    R_z( 1, 2, : ) = - sin( theta );
    R_z( 2, 1, : ) = - R_z( 1, 2, : );
    R_z( 3, 3, : ) = 1;
    
    % rotation matrices of par.n_inter random RF pulses
    
    alpha = pi .* ( rand( par.n_inter, 1 ) - 0.5 );
    phase = ( 2 * pi ) .* ( rand( par.n_inter, 1 ) - 0.5 );
    
    R_x = zeros( 3, 3, par.n_inter );
    
    R_x( 1, 1, : ) = 1;
    R_x( 2, 2, : ) = cos( alpha );
    R_x( 3, 3, : ) = R_x( 2, 2, : );
    R_x( 2, 3, : ) = - sin( alpha );
    R_x( 3, 2, : ) = - R_x( 2, 3, : );
    
    R_ph = zeros( 3, 3, par.n_inter );
    
    R_ph( 1, 1, : ) = cos( phase );
    R_ph( 2, 2, : ) = R_ph( 1, 1, : );
    R_ph( 1, 2, : ) = - sin( phase );
    R_ph( 2, 1, : ) = - R_ph( 1, 2, : );
    R_ph( 3, 3, : ) = 1;
    
    R_rf = zeros( 3, 3, par.n_inter );
    
    for i = 1 : par.n_inter
        
        R_rf( :, :, i ) = R_ph( :, :, i ) * R_x( :, :, i ) * R_ph( :, :, i ).';
        
    end
    
    % repolarization per TR
    
    m_new = zeros( 3, d );
    m_new( 3, : ) = 1 - E1;
    
    % initial state: longitudinal magnetization
    
    m = zeros( 3, 1 );
    m( 3, : ) = 1;

    % real valued magnetization (after RF pulse)
        
    m_bloch = zeros( 3, par.n_inter );
           
    %% initialize configuration model
    
    % create instance
    
    cm_exact = CoMoTk;
    
    % mandatory tissue parameters
    
    cm_exact.R1 = 1 / par.T1;
    cm_exact.R2 = 1 / par.T2;
    cm_exact.D = 0;
    
    % get default options
    
    options = cm_exact.options;
    
    % adapt options
    
    options.alloc_d = d;
    options.epsilon = 0;
    
    % set new options
    
    cm_exact.options = options;
    
    % timing
    
    time_exact = zeros( 1, par.n_inter );
    
    % real valued magnetization (after RF pulse)
    
    m_exact = zeros( 3, par.n_inter );
    
    % number of stored configurations
    
    n_conf_exact = zeros( 1, par.n_inter );
    
    % absolute deviation from Bloch simulation
    
    abs_err_exact = zeros( 1, par.n_inter );
        
    % create instance
    
    cm_rapid = CoMoTk;
    
    % mandatory tissue parameters
    
    cm_rapid.R1 = 1 / par.T1;
    cm_rapid.R2 = 1 / par.T2;
    cm_rapid.D = 0;
    
    % get default options
    
    options = cm_rapid.options;
    
    % adapt options
    
    options.alloc_d = d;
    options.epsilon = par.epsilon;
    options.rapid_meltdown = true;
    
    % set new options
    
    cm_rapid.options = options;
    
    % timing
    
    time_rapid = zeros( 1, par.n_inter );
    
    % real valued magnetization (after RF pulse)
    
    m_rapid = zeros( 3, par.n_inter );
    
    % number of stored configurations
    
    n_conf_rapid = zeros( 1, par.n_inter );
    
    % absolute deviation from Bloch simulation
    
    abs_err_rapid = zeros( 1, par.n_inter );
        
    % randomly select par.n_inter intervals from set { 1, 2 }
    
    i_tau = randi( d, par.n_inter, 1 );

    % but always start with 2 (just to simplify the code below for displaying 
    % the occupied configurations a bit)
    
    i_tau( 1 ) = 2;
    
    % count number of intervals for each dimension
    % needed to show occupied configuration orders
    
    [ ~, ~, idx_ ] = unique( i_tau );
    nd_ = accumarray( idx_, 1 )';
    
    % define matrices for occupied and stored configurations
    
    nx = - nd_( 1 ) : nd_( 1 );
    nx0 = nd_( 1 )  + 1;
    Nx = 2 .* nd_( 1 ) + 1;
    xticks_ = nx0 - 100 * floor( nd_( 1 ) ) : 100 : nx0 + 100 * floor( nd_( 1 ) );
    xticklabels_ = split( num2str( xticks_ - nx0 ) );

    ny = - nd_( 2 ) : nd_( 2 );
    ny0 = nd_( 2 )  + 1;
    Ny = 2 .* nd_( 2 ) + 1;
    yticks_ = ny0 - 100 * floor( nd_( 2 ) ) : 100 : ny0 + 100 * floor( nd_( 2 ) );
    yticklabels_ = split( num2str( ny0 - yticks_ ) );
    
    N = [ Ny, Nx ];
    
    m_ex = zeros( N );
    m_ra = zeros( N );
    
    % extract complete CM sum (like in bSSFP)
    
    sum_par = [];
    sum_par.omega = omega;
    
    %% start actual calculations
    
    for i = 1 : par.n_inter
    
        %% Bloch simulation
        
        m = R_rf( :, :, i ) * m;                                           % RF pulse
        
        m_bloch( :, i ) = m;                                               % save magnetization vector
        
        m = R_z( :, :, i_tau( i ) ) * E( :, :, i_tau( i ) ) * m + ...      % precession, relaxation
            m_new( :, i_tau( i ) );                                        % and repolarization
        
        %% Configuration model
        
        % RF parameters
        
        RF_par = [];
        RF_par.FlipAngle = alpha( i );
        RF_par.Phase = phase( i );

        % precession, relaxation and repolarization
        
        Time_par = [];
        Time_par.lambda = i_tau( i );
        Time_par.tau = tau( i_tau( i ) );
        
        % execute sequence
        
        tic;
        
        % execute RF pulse
        
        cm_exact.RF( RF_par );
        
        % calculate sum over all configurations
        
        res = cm_exact.sum( sum_par );
        
        % time interval
        
        cm_exact.time( Time_par );
        
        % store time per TR
        
        time_exact( i ) = toc;
        
        % save results
        
        m_exact( 1, i ) = real( res.xy );
        m_exact( 2, i ) = imag( res.xy );
        m_exact( 3, i ) = res.z;
        
        % absolute deviation from Bloch (complete vector)
        
        abs_err_exact( i ) = sqrt( sum( abs( m_bloch( :, i ) - m_exact( :, i ) ).^2 ) );
        
        % number of occupied (== stored) configurations
        
        n_conf_exact( i ) = cm_exact.n_conf;
        
        tic;
        
        % RF pulse
        
        cm_rapid.RF( RF_par );
        
        % calculate sum over all configurations
        
        res = cm_rapid.sum( sum_par );
        
        % time interval
        
        cm_rapid.time( Time_par );
        
        % store time per TR
        
        time_rapid( i ) = toc;
        
        % save results
        
        m_rapid( 1, i ) = real( res.xy );
        m_rapid( 2, i ) = imag( res.xy );
        m_rapid( 3, i ) = res.z;
        
        % absolute deviation from Bloch (complete vector)
        
        abs_err_rapid( i ) = sqrt( sum( abs( m_bloch( :, i ) - m_rapid( :, i ) ).^2 ) );
        
        % number of stored configurations
        
        n_conf_rapid( i ) = cm_rapid.n_conf;
        
        if ( mod( i, 10 ) == 0 || i == par.n_inter )
            
            % show occupied vs stored configurations during simulation
            
            n_occ_ex = cm_exact.n( cm_exact.b_n, : );
            i_occ_ex = n_occ_ex + [ ny0, nx0 ];
            
            n_occ_ra = cm_rapid.n( cm_rapid.b_n, : );
            i_occ_ra = n_occ_ra + [ ny0, nx0 ];
            
            co = colormap( 'parula' );
            co( 1, : ) = 0;
            
            min_ = - 20;
            
            m_occ_ex = log( abs( cm_exact.m( 1, cm_exact.b_n ) ) );
            m_occ_ex( m_occ_ex <= min_ ) = min_ + 100 * eps;
            m_ex( : ) = 0;
            m_ex( sub2ind( N,  i_occ_ex( :, 1 ), i_occ_ex( :, 2 ) ) ) = m_occ_ex - min_;
            m_ex = ceil( ( 63 / max( m_ex( : ) ) ) .* m_ex );
            m_ex( end : -1 : 1, : ) = m_ex;
            
            m_occ_ra = log( abs( cm_rapid.m( 1, cm_rapid.b_n ) ) );
            m_occ_ra( m_occ_ra <= min_ ) = min_ + 100 * eps;
            m_ra( : ) = 0;
            m_ra( sub2ind( N,  i_occ_ra( :, 1 ), i_occ_ra( :, 2 ) ) ) = m_occ_ra - min_;
            m_ra = ceil( ( 63 / max( m_ra( : ) ) ) .* m_ra );
            m_ex( end : -1 : 1, : ) = m_ex;
            
        end
            
        if ( i < par.n_inter && mod( i, 10 ) == 0 )
            
            fprintf( 1, '%d / %d\n', i, par.n_inter );
            
            subplot( 1, 2, 1 );
            imagesc( m_ex, [ 0, 64 ] );
            title( 'occupied configurations' );
            xticks( xticks_ );
            xticklabels( xticklabels_ );
            xlabel( '$n_1$', 'Interpreter', 'latex' );
            yticks( yticks_ );
            yticklabels( yticklabels_ );
            ylabel( '$n_2$', 'Interpreter', 'latex' );
            
            subplot( 1, 2, 2 );
            imagesc( m_ra, [ 0, 64 ] );
            title( 'stored configurations' );
            xticks( xticks_ );
            xticklabels( xticklabels_ );
            xlabel( '$n_1$', 'Interpreter', 'latex' );
            yticks( yticks_ );
            yticklabels( yticklabels_ );
            ylabel( '$n_2$', 'Interpreter', 'latex' );
            
            colormap( co );
            
            drawnow;
            
        elseif ( i == par.n_inter )

            % final plot for publication
            
            % some parameters
            
            rng = ( - 10 : 0 ) + par.n_inter;
            
            m_xy_bloch = m_bloch( 1, rng ) + 1i .* m_bloch( 2, rng );
            m_xy_exact = m_exact( 1, rng ) + 1i .* m_exact( 2, rng );
            m_xy_rapid = m_rapid( 1, rng  ) + 1i .* m_rapid( 2, rng );

            abs_err_exact = sqrt( sum( abs( m_exact - m_bloch ).^2, 1 ) );
            abs_err_rapid = sqrt( sum( abs( m_rapid - m_bloch ).^2, 1 ) );
            
            rng_tot = 1 : par.n_inter;
            
            % now the plots

            % size in [cm]
            
            width = 18;
            height = 12;
            
            % #1
            
            ax = subplot( 2, 3, 1 );

            % needed once
            
            set( gcf, 'Units', 'centimeters' );
            set( gcf, 'Position', [ 0, 0, width, height ] );
            set( gcf, 'Color', 'w' );
            
            delta = 0.01;
        
            %%

            plot( rng_tot, alpha, '.', rng_tot, phase, '.' ); 
            ylim( [ - pi, pi ] );
            ang_ = [ - pi, - 0.5 * pi, 0, 0.5 * pi, pi ];
            yticks( ang_ );
            yticklabels( { '-\pi', '-\pi/2', '0', '\pi/2', '\pi' } ) %, 'Interpreter', 'latex' );
            ylabel( '$\alpha_\nu, \varphi_\nu$', 'Interpreter', 'latex' );            
            xlabel( '$\nu$', 'Interpreter', 'latex' );
            xlim( [ 0, par.n_inter ] );
            legend( 'flip angle', 'phase', 'Interpreter', 'latex', 'Location', 'south' );

            title( 'Random Sequence', 'Interpreter', 'latex' );
                
            ti = ax.TightInset;
            left = ti(1) + delta;
            bottom = ti(2) + delta + 0.5;
            ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            ax.Position = [left bottom ax_width ax_height ];
            
            %%
            
            % #2

            ax = subplot( 2, 3, 2 );
            
            plot( rng, abs( m_xy_bloch ), '+k', rng, abs( m_xy_exact ), 'b', rng, abs( m_xy_rapid ), 'r' );
            xlabel( '$\nu$', 'Interpreter', 'latex' );
            ylabel( '$\left|m_{xy}\right|$', 'Interpreter', 'latex' );
            legend( 'Bloch', 'CM (exact)', 'CM (rapid)', 'Interpreter', 'latex', 'Location', 'best' );
            title( 'Transverse Magnetization', 'Interpreter', 'latex' );
            
            ti = ax.TightInset;
            left = 0.3333 + ti(1) + delta;
            bottom = ti(2) + delta + 0.5;
            ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            ax.Position = [left bottom ax_width ax_height ];

            % #3
            
            ax = subplot( 2, 3, 3 );

            semilogy( rng_tot, abs_err_exact, rng_tot, abs_err_rapid );
            xlim( [ 0 par.n_inter ] );
            xlabel( '$\nu$', 'Interpreter', 'latex' );
            ylabel( '$\left|\bf{\Delta m}\right|$', 'Interpreter', 'latex' );            
            legend( 'exact', 'rapid', 'Interpreter', 'latex', 'Location', 'east' );
            title( 'Deviation from Bloch', 'Interpreter', 'latex' );
            
            ti = ax.TightInset;
            left = 0.6667 + ti(1) + delta;
            bottom = ti(2) + delta + 0.5;
            ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            ax.Position = [left bottom ax_width ax_height ];

            % #4
            
            ax = subplot( 2, 3, 4 );

            semilogy( rng_tot, n_conf_exact, rng_tot, n_conf_rapid );
            xlim( [ 0 par.n_inter ] );
            xlabel( '$\nu$', 'Interpreter', 'latex' );
            ylabel( 'number', 'Interpreter', 'latex' );
            legend( 'occupied', 'stored', 'Interpreter', 'latex', 'Location', 'south' );
            title( 'Configurations', 'Interpreter', 'latex' );

            ti = ax.TightInset;
            left = ti(1) + delta;
            bottom = ti(2) + delta;
            ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            ax.Position = [left bottom ax_width ax_height ];

            % #5

            ax = subplot( 2, 3, 5 );

            imagesc( m_ra, [ 0, 64 ] );
            xticks( xticks_ );
            xticklabels( xticklabels_ );
            xlabel( '$n_1$', 'Interpreter', 'latex' );
            yticks( yticks_ );
            yticklabels( yticklabels_ );
            ylabel( '$n_2$', 'Interpreter', 'latex' );
            title( 'Stored Configurations', 'Interpreter', 'latex' );
            
            ti = ax.TightInset;
            left = 0.3333 + ti(1) + delta;
            bottom = ti(2) + delta;
            ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            ax.Position = [left bottom ax_width ax_height ];

            % #6
                        
            ax = subplot( 2, 3, 6 );
            
            plot( rng_tot, time_exact, rng_tot, time_rapid );
            xlim( [ 0 par.n_inter ] );
            xlabel( '$\nu$', 'Interpreter', 'latex' );
            ylabel( '[ms]', 'Interpreter', 'latex' );
            legend( 'exact', 'rapid', 'Interpreter', 'latex', 'Location', 'northwest' );
            title( 'CPU Time / Interval', 'Interpreter', 'latex' );
            
            ti = ax.TightInset;
            left = 0.6667 + ti(1) + delta;
            bottom = ti(2) + delta;
            ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            ax.Position = [left bottom ax_width ax_height ];

            colormap( co );
            
        end
        
    end
    
    
end
