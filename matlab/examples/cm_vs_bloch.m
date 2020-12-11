%% Comparison CoMoTk vs Bloch simulation for a random RF train
% 
% - no gradients, i.e. all magnetization pathways contribute
% - flip angle and phases randomized for each RF pulse
% - variable time intervals
% - single-peak tissue with freely adjustable off-resonance
% - no susceptibility effects (i.e. microscopic solution)

par = [];
opt = [];
str = [];

par.T1 = 100;
par.T2 = 10;
par.tau_min = 0.8;
par.tau_max = 1.2;
par.n_tau = 10;
par.omega = 1;
par.n_inter = 100;

opt.T1 = [];
opt.T2 = [];
opt.tau_min = [];
opt.tau_max = [];
opt.n_tau = [];
opt.omega = [];
opt.n_inter = [];

str.T1 = '[ms] T1';
str.T2 = '[ms] T2';
str.tau_min = '[ms] minimum duration between RF pulses.';
str.tau_max = '[ms] maximum duration between RF pulses.';
str.n_tau = '[ms] allowed RF pulse spacing: linspace(tau_min,tau_max,n_tau)';
str.omega = '[rad/ms] dephasing angle due to off-resonance (randomly chosen, if empty)';
str.n_inter = 'number of RF pulses (separated by non-equidistant time intervals)';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
       
    % allowed times
    
    n_tau = par.n_tau;
    tau = linspace( par.tau_min, par.tau_max, n_tau );
    d_tau = tau( 2 ) - tau( 1 );
    
    % off-resonance frequency
    
    if ( isempty( par.omega ) )
    
        % random setting, if not specified
        
        omega = 2 * pi * rand;
        
    else
        
        % otherwise set by the user
        
        omega = par.omega;
        
    end
    
    %% initialize Bloch simulations
    
    % dephasing angle due to off-resonance, depending on the n_tau durations tau(:)
    
    theta = - omega .* tau;
    
    % n_tau relaxation and precession precession matrices
    
    E = zeros( 3, 3, n_tau );
    
    E1 = exp( - tau ./ par.T1 );
    E2 = exp( - tau ./ par.T2 );
    
    E( 1, 1, : ) = E2;
    E( 2, 2, : ) = E2;
    E( 3, 3, : ) = E1;
    
    R_z = zeros( 3, 3, n_tau );
    
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
    
    m_new = zeros( 3, n_tau );
    m_new( 3, : ) = 1 - E1;
    
    % initial state: longitudinal magnetization
    
    m = zeros( 3, 1 );
    m( 3, : ) = 1;

    % real valued magnetization (after RF pulse)
        
    m_bloch = zeros( 3, par.n_inter );
           
    %% initialize configuration model
    
    % create instance
    
    cm = CoMoTk;
    
    % mandatory tissue parameters
    
    cm.R1 = 1 / par.T1;
    cm.R2 = 1 / par.T2;
    cm.D = 0;
      
    % support in (p,tau) space
    % no gradients are applied, therefore we only need to set the time
    % dimension
    
    cm.d_tau = d_tau;
    cm.n_tau = par.n_inter * par.tau_max / d_tau; % should be a safe choice
    
    %% other settings
    % timing
    
    time_cm = zeros( 1, par.n_inter );
    
    % real valued magnetization (after RF pulse)
    
    m_cm = zeros( 3, par.n_inter );
    
    % total and fractional number of stored configurations
    
    occ = zeros( 1, par.n_inter );
    fra = zeros( 1, par.n_inter );
    
    % absolute deviation from Bloch simulation
    
    abs_err = zeros( 1, par.n_inter );
        
    % randomly select par.n_inter intervals from set { 1, ..., n_tau }
    
    i_tau = randi( n_tau, par.n_inter, 1 );
    
    % extract complete CM sum (like in bSSFP)
    
    sum_par = [];

    if ( omega ~= 0 )
        
        sum_par.omega = omega;
    
    end
        
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
        Time_par.tau = tau( i_tau( i ) );
        
        % execute sequence
        
        tic;
        
        % execute RF pulse
        
        cm.RF( RF_par );
        
        % calculate sum over all configurations
        
        res = cm.sum( sum_par );
        
        % time interval
        
        cm.time( Time_par );
        
        % store time per TR
        
        time_cm( i ) = toc;
        
        % save results
        
        m_cm( 1, i ) = real( res.xy );
        m_cm( 2, i ) = imag( res.xy );
        m_cm( 3, i ) = res.z;
        
        % absolute deviation from Bloch (complete vector)
        
        abs_err( i ) = sqrt( sum( abs( m_bloch( :, i ) - m_cm( :, i ) ).^2 ) );
        
        % number of occupied (== stored) configurations in (p,tau) space
        
        [ o_, f_ ] = cm.n_conf;
        occ( i ) = o_;
        fra( i ) = f_( 4 ); % corresponding to the tau direction
        
        if ( mod( i, 10 ) == 0 || i == par.n_inter )
            
            % show occupied vs stored configurations during simulation
            
            %             n_occ_ex = cm.n( cm.b_n, : );
            %             i_occ_ex = n_occ_ex + [ ny0, nx0 ];
            %
            %             n_occ_ra = cm_rapid.n( cm_rapid.b_n, : );
            %             i_occ_ra = n_occ_ra + [ ny0, nx0 ];
            %
            %             co = colormap( 'parula' );
            %             co( 1, : ) = 0;
            %
            %             min_ = - 20;
            %
            %             m_occ_ex = log( abs( cm.m( 1, cm.b_n ) ) );
            %             m_occ_ex( m_occ_ex <= min_ ) = min_ + 100 * eps;
            %             m_ex( : ) = 0;
            %             m_ex( sub2ind( N,  i_occ_ex( :, 1 ), i_occ_ex( :, 2 ) ) ) = m_occ_ex - min_;
            %             m_ex = ceil( ( 63 / max( m_ex( : ) ) ) .* m_ex );
            %             m_ex( end : -1 : 1, : ) = m_ex;
            %
            %             m_occ_ra = log( abs( cm_rapid.m( 1, cm_rapid.b_n ) ) );
            %             m_occ_ra( m_occ_ra <= min_ ) = min_ + 100 * eps;
            %             m_ra( : ) = 0;
            %             m_ra( sub2ind( N,  i_occ_ra( :, 1 ), i_occ_ra( :, 2 ) ) ) = m_occ_ra - min_;
            %             m_ra = ceil( ( 63 / max( m_ra( : ) ) ) .* m_ra );
            %             m_ex( end : -1 : 1, : ) = m_ex;
            
        end
            
        if ( i < par.n_inter && mod( i, 10 ) == 0 )
            
            fprintf( 1, '%d / %d\n', i, par.n_inter );
            
            %             subplot( 1, 2, 1 );
            %             imagesc( m_ex, [ 0, 64 ] );
            %             title( 'occupied configurations' );
            %             xticks( xticks_ );
            %             xticklabels( xticklabels_ );
            %             xlabel( '$n_1$', 'Interpreter', 'latex' );
            %             yticks( yticks_ );
            %             yticklabels( yticklabels_ );
            %             ylabel( '$n_2$', 'Interpreter', 'latex' );
            %
            %             subplot( 1, 2, 2 );
            %             imagesc( m_ra, [ 0, 64 ] );
            %             title( 'stored configurations' );
            %             xticks( xticks_ );
            %             xticklabels( xticklabels_ );
            %             xlabel( '$n_1$', 'Interpreter', 'latex' );
            %             yticks( yticks_ );
            %             yticklabels( yticklabels_ );
            %             ylabel( '$n_2$', 'Interpreter', 'latex' );
            %
            %             colormap( co );
            %
            %             drawnow;
            
        elseif ( i == par.n_inter )

            % final plot for publication
            
            % some parameters
            
            %             rng = ( - 10 : 0 ) + par.n_inter;
            %
            %             m_xy_bloch = m_bloch( 1, rng ) + 1i .* m_bloch( 2, rng );
            %             m_xy_exact = m_cm( 1, rng ) + 1i .* m_cm( 2, rng );
            %             m_xy_rapid = m_rapid( 1, rng  ) + 1i .* m_rapid( 2, rng );
            %
            %             abs_err = sqrt( sum( abs( m_cm - m_bloch ).^2, 1 ) );
            %             abs_err_rapid = sqrt( sum( abs( m_rapid - m_bloch ).^2, 1 ) );
            %
            %             rng_tot = 1 : par.n_inter;
            %
            %             % now the plots
            %
            %             % size in [cm]
            %
            %             width = 18;
            %             height = 12;
            %
            %             % #1
            %
            %             ax = subplot( 2, 3, 1 );
            %
            %             % needed once
            %
            %             set( gcf, 'Units', 'centimeters' );
            %             set( gcf, 'Position', [ 0, 0, width, height ] );
            %             set( gcf, 'Color', 'w' );
            %
            %             delta = 0.01;
            %
            %             %%
            %
            %             plot( rng_tot, alpha, '.', rng_tot, phase, '.' );
            %             ylim( [ - pi, pi ] );
            %             ang_ = [ - pi, - 0.5 * pi, 0, 0.5 * pi, pi ];
            %             yticks( ang_ );
            %             yticklabels( { '-\pi', '-\pi/2', '0', '\pi/2', '\pi' } ) %, 'Interpreter', 'latex' );
            %             ylabel( '$\alpha_\nu, \varphi_\nu$', 'Interpreter', 'latex' );
            %             xlabel( '$\nu$', 'Interpreter', 'latex' );
            %             xlim( [ 0, par.n_inter ] );
            %             legend( 'flip angle', 'phase', 'Interpreter', 'latex', 'Location', 'south' );
            %
            %             title( 'Random Sequence', 'Interpreter', 'latex' );
            %
            %             ti = ax.TightInset;
            %             left = ti(1) + delta;
            %             bottom = ti(2) + delta + 0.5;
            %             ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            %             ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            %             ax.Position = [left bottom ax_width ax_height ];
            %
            %             %%
            %
            %             % #2
            %
            %             ax = subplot( 2, 3, 2 );
            %
            %             plot( rng, abs( m_xy_bloch ), '+k', rng, abs( m_xy_exact ), 'b', rng, abs( m_xy_rapid ), 'r' );
            %             xlabel( '$\nu$', 'Interpreter', 'latex' );
            %             ylabel( '$\left|m_{xy}\right|$', 'Interpreter', 'latex' );
            %             legend( 'Bloch', 'CM (exact)', 'CM (rapid)', 'Interpreter', 'latex', 'Location', 'best' );
            %             title( 'Transverse Magnetization', 'Interpreter', 'latex' );
            %
            %             ti = ax.TightInset;
            %             left = 0.3333 + ti(1) + delta;
            %             bottom = ti(2) + delta + 0.5;
            %             ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            %             ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            %             ax.Position = [left bottom ax_width ax_height ];
            %
            %             % #3
            %
            %             ax = subplot( 2, 3, 3 );
            %
            %             semilogy( rng_tot, abs_err, rng_tot, abs_err_rapid );
            %             xlim( [ 0 par.n_inter ] );
            %             xlabel( '$\nu$', 'Interpreter', 'latex' );
            %             ylabel( '$\left|\bf{\Delta m}\right|$', 'Interpreter', 'latex' );
            %             legend( 'exact', 'rapid', 'Interpreter', 'latex', 'Location', 'east' );
            %             title( 'Deviation from Bloch', 'Interpreter', 'latex' );
            %
            %             ti = ax.TightInset;
            %             left = 0.6667 + ti(1) + delta;
            %             bottom = ti(2) + delta + 0.5;
            %             ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            %             ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            %             ax.Position = [left bottom ax_width ax_height ];
            %
            %             % #4
            %
            %             ax = subplot( 2, 3, 4 );
            %
            %             semilogy( rng_tot, occ, rng_tot, n_conf_rapid );
            %             xlim( [ 0 par.n_inter ] );
            %             xlabel( '$\nu$', 'Interpreter', 'latex' );
            %             ylabel( 'number', 'Interpreter', 'latex' );
            %             legend( 'occupied', 'stored', 'Interpreter', 'latex', 'Location', 'south' );
            %             title( 'Configurations', 'Interpreter', 'latex' );
            %
            %             ti = ax.TightInset;
            %             left = ti(1) + delta;
            %             bottom = ti(2) + delta;
            %             ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            %             ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            %             ax.Position = [left bottom ax_width ax_height ];
            %
            %             % #5
            %
            %             ax = subplot( 2, 3, 5 );
            %
            %             imagesc( m_ra, [ 0, 64 ] );
            %             xticks( xticks_ );
            %             xticklabels( xticklabels_ );
            %             xlabel( '$n_1$', 'Interpreter', 'latex' );
            %             yticks( yticks_ );
            %             yticklabels( yticklabels_ );
            %             ylabel( '$n_2$', 'Interpreter', 'latex' );
            %             title( 'Stored Configurations', 'Interpreter', 'latex' );
            %
            %             ti = ax.TightInset;
            %             left = 0.3333 + ti(1) + delta;
            %             bottom = ti(2) + delta;
            %             ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            %             ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            %             ax.Position = [left bottom ax_width ax_height ];
            %
            %             % #6
            %
            %             ax = subplot( 2, 3, 6 );
            %
            %             plot( rng_tot, time_cm, rng_tot, time_rapid );
            %             xlim( [ 0 par.n_inter ] );
            %             xlabel( '$\nu$', 'Interpreter', 'latex' );
            %             ylabel( '[ms]', 'Interpreter', 'latex' );
            %             legend( 'exact', 'rapid', 'Interpreter', 'latex', 'Location', 'northwest' );
            %             title( 'CPU Time / Interval', 'Interpreter', 'latex' );
            %
            %             ti = ax.TightInset;
            %             left = 0.6667 + ti(1) + delta;
            %             bottom = ti(2) + delta;
            %             ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
            %             ax_height = 0.5 - ti(2) - ti(4) - 2 * delta;
            %             ax.Position = [left bottom ax_width ax_height ];
            %
            %             colormap( co );
            
        end
        
    end
    
    
end
