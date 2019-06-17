%% Comparison CoMoTk vs Bloch simulation for a random RF train

par = [];
opt = [];
str = [];

par.T1 = 100;
par.T2 = 10;
par.D = 0;
par.n_inter = 50;
par.d = 5;
par.n_RF = 5;
par.epsilon = 1e-4;
par.alloc_n = 1000;
par.alloc_d = 2;
par.rapid_meltdown = 'True';
par.verbose = 'False';

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.n_inter = [];
opt.d = [];
opt.n_RF = [];
opt.epsilon = [];
opt.alloc_n = [];
opt.alloc_d = [];
opt.rapid_meltdown = { 'True', 'False' };
opt.verbose = { 'True', 'False' };

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.n_inter = 'number of RF pulses (separated by non-equidistant time intervals)';
str.d = 'dimension of configuration model (= number of different time intervals)';
str.n_RF = 'number of different RF pulses (< n_inter)';
str.epsilon = 'discard configuration vectors with L2 norm smaller than this';
str.alloc_n = 'initially allocated configurations';
str.alloc_d = 'initially allocated dimensions';
str.rapid_meltdown = 'faster, but more memory';
str.verbose = 'provide some informal output';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
        
    % The RF pulses are selected from a subset of par.n_RF < par.n_inter random pulses.
    
    alpha = 0.5 .* pi .* randn( par.n_RF, 1 );
    phase = pi .* randn( par.n_RF, 1 );
    
    % The time intervals are chosen from a set of par.d < par.n_inter random durations (between 0.5 and 1.5).
    
    tau = 0.5 + rand( par.d, 1 );
    
    % off-resonance frequency and associated accumulated phase in the par.d intervals
    
    omega = pi;
    theta = omega .* tau;
    
    %% initialize Bloch simulation
    
    % par.d relaxation and precession precession matrices
    
    E = zeros( 3, 3, par.d );
    
    E1 = exp( - tau / par.T1 );
    E2 = exp( - tau / par.T2 );
    
    E( 1, 1, : ) = E2;
    E( 2, 2, : ) = E2;
    E( 3, 3, : ) = E1;
    
    R_z = zeros( 3, 3, par.d );
    
    R_z( 1, 1, : ) = cos( theta );
    R_z( 2, 2, : ) = R_z( 1, 1, : );
    R_z( 1, 2, : ) = - sin( theta );
    R_z( 2, 1, : ) = - R_z( 1, 2, : );
    R_z( 3, 3, : ) = 1;
    
    % par.n_RF rotation matrices
    
    R_x = zeros( 3, 3, par.n_RF );
    
    R_x( 1, 1, : ) = 1;
    R_x( 2, 2, : ) = cos( alpha );
    R_x( 3, 3, : ) = R_x( 2, 2, : );
    R_x( 2, 3, : ) = - sin( alpha );
    R_x( 3, 2, : ) = - R_x( 2, 3, : );
    
    R_ph = zeros( 3, 3, par.n_RF );
    
    R_ph( 1, 1, : ) = cos( phase );
    R_ph( 2, 2, : ) = R_ph( 1, 1, : );
    R_ph( 1, 2, : ) = - sin( phase );
    R_ph( 2, 1, : ) = - R_ph( 1, 2, : );
    R_ph( 3, 3, : ) = 1;
    
    R_rf = zeros( 3, 3, par.n_RF );
    
    for i = 1 : par.n_RF
        
        R_rf( :, :, i ) = R_ph( :, :, i ) * R_x( :, :, i ) * R_ph( :, :, i ).';
        
    end
    
    % repolarization per TR
    
    m_new = zeros( 3, par.d );
    m_new( 3, : ) = 1 - E1;
    
    % initial state: longitudinal magnetization
    
    m = zeros( 3, 1 );
    m( 3, : ) = 1;

    % set meltdown algorithm
    
    if ( isequal( par.rapid_meltdown, 'True' ) )
        
        rapid_meltdown = true;
        
    else
        
        rapid_meltdown = false;
        
    end
    
    % set verbosity
    
    if ( isequal( par.verbose, 'True' ) )
        
        verbose = true;
        
    else
        
        verbose = false;
        
    end
        
    
    %% initialize configuration model
    
    cm = CoMoTk;
    
    % mandatory tissue parameters
    
    cm.R1 = 1 / par.T1;
    cm.R2 = 1 / par.T2;
    cm.D = par.D;
    
    % get default options
    
    options = cm.options;
    
    % adapt options
    
    options.alloc_n = par.alloc_n;
    options.alloc_d = par.alloc_d;
    options.epsilon = par.epsilon;
    options.rapid_meltdown = rapid_meltdown;
    options.verbose = verbose;
    
    % set new options
    
    cm.options = options;
    
    % initialize everything
    
    cm.init_configuration ( [ 0; 0; 1 ] );
    
    %% allocate space for results
    
    m_bloch_p = zeros( 3, par.n_inter );     % after RF pulse
    m_bloch_m = zeros( 3, par.n_inter );    % before RF pulse
    
    m_iso_p = zeros( 3, par.n_inter );     % after RF pulse
    m_iso_m = zeros( 3, par.n_inter );    % before RF pulse
    
    % select RF pulses and times randomly
    
    i_rf = randi( par.n_RF, par.n_inter, 1 );
    i_tau = randi( par.d, par.n_inter, 1 );
    
    for i = 1 : par.n_inter
    
        tic;
    
        fprintf( 1, '%d / %d\n', i, par.n_inter );
        
        % Bloch simulation
        
        m = R_rf( :, :, i_rf( i ) ) * m;                                      % RF pulse
        
        m_bloch_p( :, i ) = m;                                                % save magnetization vector
        
        m = R_z( :, :, i_tau( i ) ) * E( :, :, i_tau( i ) ) * m + ...         % precession, relaxation
            m_new( :, i_tau( i ) );                                           % and repolarization
        
        m_bloch_m( :, i ) = m;
        
        % Configuration model
        
        param = [];
        param.FlipAngle = alpha( i_rf( i ) );
        param.Phase = phase( i_rf( i ) );
        
        cm.RF( param );
        
        %   cm.check_state();
        
        param = [];
        param.omega = omega;
        
        res = cm.sum( param );                                 % get sum
        
        m_iso_p( 1, i ) = real( res.xy );                                     % store as real vector
        m_iso_p( 2, i ) = imag( res.xy );
        m_iso_p( 3, i ) = res.z;
        
        % precession, relaxation and repolarization
        
        param = [];
        param.mu = i_tau( i );
        param.tau = tau( i_tau( i ) );
        
        cm.time( param );
        
        %    cm.check_state();
        
        param = [];
        param.omega = omega;
        
        res = cm.sum( param );                                 % get sum
        
        m_iso_m( 1, i ) = real( res.xy );                                     % store as real vector
        m_iso_m( 2, i ) = imag( res.xy );
        m_iso_m( 3, i ) = res.z;
        
        fprintf( 1, 'After RF  : Deviation = %e\n', ...
            sqrt( sum( abs( m_bloch_p( :, i ) - m_iso_p( :, i ) ).^2 ) ) );
        fprintf( 1, 'After Time: Deviation = %e\n', ...
            sqrt( sum( abs( m_bloch_m( :, i ) - m_iso_m( :, i ) ).^2 ) ) );

        toc;
    
    end
    
end
