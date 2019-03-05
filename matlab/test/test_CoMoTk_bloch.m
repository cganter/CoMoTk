%% Comparison CoMoTk vs Bloch simulation for a random RF train

% Tissue

T1 = 1000;
T2 = 100;
D = 0;

% - We consider a sequence of n instantaneous RF pulses, separated by non-equidistant time intervals.

n = 100;

% The RF pulses are selected from a subset of n_rf < n random pulses.

n_rf = 5;
alpha = 0.5 .* pi .* randn( n_rf, 1 );
phase = pi .* randn( n_rf, 1 );

% The time intervals are chosen from a set of n_tau < n random durations (between 0.5 and 1.5).

n_tau = 3;
tau = 0.5 + rand( n_tau, 1 );

% off-resonance frequency and associated accumulated phase in the n_tau intervals

omega = pi;
theta = omega .* tau;

%% initialize Bloch simulation

% n_tau relaxation and precession precession matrices

E = zeros( 3, 3, n_tau );

E1 = exp( - tau / T1 );
E2 = exp( - tau / T2 );

E( 1, 1, : ) = E2;
E( 2, 2, : ) = E2;
E( 3, 3, : ) = E1;

R_z = zeros( 3, 3, n_tau );

R_z( 1, 1, : ) = cos( theta );
R_z( 2, 2, : ) = R_z( 1, 1, : );
R_z( 1, 2, : ) = - sin( theta );
R_z( 2, 1, : ) = - R_z( 1, 2, : );
R_z( 3, 3, : ) = 1;

% n_rf rotation matrices

R_x = zeros( 3, 3, n_rf );

R_x( 1, 1, : ) = 1;
R_x( 2, 2, : ) = cos( alpha );
R_x( 3, 3, : ) = R_x( 2, 2, : );
R_x( 2, 3, : ) = - sin( alpha );
R_x( 3, 2, : ) = - R_x( 2, 3, : );

R_ph = zeros( 3, 3, n_rf );

R_ph( 1, 1, : ) = cos( phase );
R_ph( 2, 2, : ) = R_ph( 1, 1, : );
R_ph( 1, 2, : ) = - sin( phase );
R_ph( 2, 1, : ) = - R_ph( 1, 2, : );
R_ph( 3, 3, : ) = 1;

R_rf = zeros( 3, 3, n_rf );

for i = 1 : n_rf
	
    R_rf( :, :, i ) = R_ph( :, :, i ) * R_x( :, :, i ) * R_ph( :, :, i ).';

end

% repolarization per TR

m_new = zeros( 3, n_tau );
m_new( 3, : ) = 1 - E1;

% initial state: longitudinal magnetization

m = zeros( 3, 1 );
m( 3, : ) = 1;

%% initialize configuration model

cm = CoMoTk;

% mandatory tissue parameters

cm.R1 = 1 / T1;
cm.R2 = 1 / T2;
cm.D = D;

% get default options

options = cm.options;

% we want to see reallocations, removal of configurations and some feedback

options.alloc_n = 5000;
options.alloc_d = 2;
options.epsilon = 0;
options.rapid_meltdown = true; % for improved speed (default == true)
options.verbose = true;
options.debug = false;

% set new options

cm.options = options;

% initialize everything

cm.init_configuration ( [ 0; 0; 1 ] );

%% allocate space for results

m_bloch_p = zeros( 3, n );     % after RF pulse
m_bloch_m = zeros( 3, n );    % before RF pulse

m_iso_p = zeros( 3, n );     % after RF pulse
m_iso_m = zeros( 3, n );    % before RF pulse

% select RF pulses and times randomly

i_rf = ceil( n_rf .* rand( n, 1 ) );
i_tau = ceil( n_tau .* rand( n, 1 ) );

tic;

for i = 1 : n

    fprintf( 1, '%d / %d\n', i, n );
    
    % Bloch simulation
    
    m = R_rf( :, :, i_rf( i ) ) * m;                                      % RF pulse
    
    m_bloch_p( :, i ) = m;                                                % save magnetization vector
    
    m = R_z( :, :, i_tau( i ) ) * E( :, :, i_tau( i ) ) * m + ...         % precession, relaxation
        m_new( :, i_tau( i ) );                                           % and repolarization
    
    m_bloch_m( :, i ) = m;
    
    % Configuration model
    
    cm.RF( alpha( i_rf( i ) ), phase( i_rf( i ) ) );                      % RF pulse
    
 %   cm.check_state();
    
    iso = cm.isochromat( omega, [], [] );                                 % get isochromat

    m_iso_p( 1, i ) = real( iso.xy );                                     % store as real vector
    m_iso_p( 2, i ) = imag( iso.xy );
    m_iso_p( 3, i ) = iso.z;
    
    % precession, relaxation and repolarization
 
    opt_time.tau = tau( i_tau( i ) );
    
    if ( ~cm.time( i_tau( i ), opt_time ) )
        
        break;
        
    end

%    cm.check_state();
        
    iso = cm.isochromat( omega, [], [] );                                 % get isochromat

    m_iso_m( 1, i ) = real( iso.xy );                                     % store as real vector
    m_iso_m( 2, i ) = imag( iso.xy );
    m_iso_m( 3, i ) = iso.z;

    fprintf( 1, 'After RF  : Deviation = %e\n', ...
             sqrt( sum( abs( m_bloch_p( :, i ) - m_iso_p( :, i ) ).^2 ) ) );
    fprintf( 1, 'After Time: Deviation = %e\n', ...
             sqrt( sum( abs( m_bloch_m( :, i ) - m_iso_m( :, i ) ).^2 ) ) );
        
end

toc;