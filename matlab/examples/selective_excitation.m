%% Simple selective excitation
% This script shows how real RF pulses can be implemented.
% The example shows a simple amplitude modulated pulse
% The generalization to phase modulated and/or adiabatic pulses works along the same lines.
% Depending on the parameters, relaxation (e.g. for T2 < t_rf) and diffusion (large ADC) 
% effects may become relevant 

par = [];
opt = [];
str = [];

par.T1 = Inf;
par.T2 = 10;
par.D = 0;
par.fa = 50;
par.sl_th = 1;
par.t_rf = 1;
par.supp_rf = 50;
par.qual_rf = 3;
par.filt_rf = 'None';
par.verbose = 'False';

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.fa = [];
opt.sl_th = [];
opt.t_rf = [];
opt.supp_rf = [];
opt.qual_rf = [];
opt.filt_rf = { 'Hamming', 'None' };
opt.verbose = { 'True', 'False' };

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.fa = '[deg] flip angle';
str.sl_th = '[mm] slice thickness';
str.t_rf = '[ms] RF pulse duration';
str.supp_rf = 'number of RF support points';
str.qual_rf = 'number of SINC pulse zero crossings (on each side)';
str.filt_rf = 'filter RF pulse to reduce wiggles';
str.verbose = 'provide some informal output';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
        
    % convert to units, as expected by CoMoTk
    
    alpha_rad = par.fa * pi / 180;
    sl_th_um = 1000 * par.sl_th;
    
    %% now we set the pulse profile
    % as described in the Handbook of MRI for a SINC pulse +/- Hamming filter
    
    % duration of small intervals
    
    tau_rf = par.t_rf / par.supp_rf;
        
    % support points and flip angles of instantaneous RF pulses (\propto B1 profile)
    
    t = par.qual_rf .* linspace( -1, 1, par.supp_rf + 1 );

    if ( isequal( par.filt_rf, 'Hamming' ) )
        
        alpha = ( 0.54 + 0.46 .* cos( pi .* t ./ par.qual_rf ) ) .* sinc( t );
        
    else
        
        alpha = sinc( t );
        
    end
    
    % normalization: at the center of the slice profile, we need the desired flip angle
    
    alpha = alpha .* ( alpha_rad / sum( alpha ) );
    
    % constant phase (purely amplitude modulated RF pulse)
    
    phase_rad = pi;
    
    %% where we want to see the slice profile
    % no off-resonance (nonzero value shifts slice)
    
    omega = 0;
    
    % locations (slice normal along z-direction)
    
    n_sl = 501;
    x = zeros( 3, n_sl );
    x( 3, : ) = linspace( - sl_th_um, sl_th_um, n_sl );
    
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
    
    options.alloc_n = 1000;
    options.alloc_d = 2;     % == slice selection and rephasing gradient
    options.epsilon = 0;     % best accuracy
    options.verbose = verbose;
    
    % set new options
    
    cm.options = options;
    
    % start with longitudinal magnetization
    
    cm.init_configuration ( [ 0; 0; 1 ] );
    
    %% calculate slice selection and rephasing gradient moment
    
    % RF pulse bandwidth
    
    f = 2 * par.qual_rf / par.t_rf;
    
    % total moment of slice selection gradient
    
    p_sl = 2 * pi * f * par.t_rf / sl_th_um;
    
    % split into par.supp_rf intervals
    
    p_rf = [ 0; 0; p_sl / par.supp_rf ];
    
    % assign unique handles for the two time intervals
    
    lambda_rf = 1;

    %% calculate SLR profile
    
    % Eq. (ll)
    
    C = cos( 0.5 .* alpha );
    S = 1i .* sin( 0.5 .* alpha );
    z05 = exp( 0.5 .* 1i .* p_rf( 3 ) .* x( 3, : ) );
    z = z05 .* z05;
    
    n = par.supp_rf + 1;
    
    A = zeros( n, n_sl );
    B = zeros( n, n_sl );
    
    A( 1, : ) = C( 1 );
    B( 1, : ) = S( 1 );
    
    for j = 2 : n
     
        A( j, : ) = C( j ) .* A( j - 1, : ) - conj( S( j ) .* z ) .* B( j - 1, : );
        B( j, : ) = S( j ) .* A( j - 1, : ) + C( j ) .* conj( z ) .* B( j - 1, : );
        
    end
    
    al_n = conj( z05.^n ) .* A( n, : );
    be_n = conj( z05.^n ) .* B( n, : );
    
    m_xy_SLR = 2 .* conj( al_n ) .* be_n;
    m_z_SLR = real( al_n .* conj( al_n ) - be_n .* conj( be_n ) );
    
    %% execute the RF pulse
    
    for i = 1 : par.supp_rf
        
        if ( i == 1 )
            
            % we start with the first small RF pulse
            
            param = [];
            param.FlipAngle = alpha( 1 );
            param.Phase = phase_rad;
            
            cm.RF( param );
            
            % in the first time interval, we specify the parameters
            
            param = [];
            param.lambda = lambda_rf;
            param.tau = tau_rf;
            param.p = p_rf;
            
            cm.time( param );
            
        else
            
            % in subsequent calls, we only need the handle
            
            param = [];
            param.lambda = lambda_rf;
            
            cm.time( param );
            
        end
        
        % and the rest of the small pulses
        
        param = [];
        param.FlipAngle = alpha( i + 1 );
        param.Phase = phase_rad;
        
        cm.RF( param );
        
    end
    
    %% collect the slice profile at locations x
    
    m_xy_CM = zeros( 1, n_sl );
    m_z_CM = zeros( 1, n_sl );
    
    for i = 1 : n_sl
        
        % we sum over all configurations (third argument == [])
        
        param = [];
        param.omega = omega;
        param.x = x( :, i );
        
        res = cm.sum( param );
 
        m_xy_CM( i ) = res.xy;
        m_z_CM( i ) = real( res.z );
        
    end
    
    %% look at the results
    
    loc = x( 3, : ) .* 0.001;
    t = par.t_rf .* linspace( -0.5, 0.5, par.supp_rf + 1 );
    
    ax = subplot( 1, 3, 1);
    stem( t, alpha ./ max( alpha ), 'filled' );
    xlim( [ - 0.5 * par.t_rf 0.5 * par.t_rf ] );
    ylim( [ min( alpha ./ max( alpha ) ) - 0.05, 1.05 ] );
    xlabel( '$t$ [ms]', 'Interpreter', 'latex' );
    ylabel( '$\alpha_\nu/\alpha_{nom}$', 'Interpreter', 'latex' );
    title( 'RF Pulse', 'Interpreter', 'latex' );

    width = 18;
    height = 6;
    
    set( gcf, 'Units', 'centimeters' );
    set( gcf, 'Position', [ 0, 0, width, height ] );
    set( gcf, 'Color', 'w' );

    delta = 0.01;
    
    ti = ax.TightInset;
    left = ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    ax = subplot( 1, 3, 2 );
    plot( loc, abs( m_xy_CM ), loc, abs( m_xy_SLR ) );
    legend( 'CM', 'SLR', 'Interpreter', 'latex', 'Location', 'south' );
    xlabel( 'position [mm]', 'Interpreter', 'latex' );
    ylabel( '$\left|m_{xy}\right|$', 'Interpreter', 'latex' );
    title( 'Transverse Part', 'Interpreter', 'latex');

    ti = ax.TightInset;
    left = 0.3333 + ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    ax = subplot( 1, 3, 3 );
    plot( loc, m_z_CM, loc, m_z_SLR );
    legend( 'CM', 'SLR', 'Interpreter', 'latex', 'Location', 'north' );
    xlabel( 'position [mm]', 'Interpreter', 'latex' );
    ylabel( '$m_z$', 'Interpreter', 'latex' );
    title( 'Longitudinal Part', 'Interpreter', 'latex' );
    
    ti = ax.TightInset;
    left = 0.6667 + ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.3333 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    fprintf( 1, 'median( abs( CM / SLR - 1 ) )\n' );
    fprintf( 1, 'xy: %e\n', median( abs( m_xy_CM / m_xy_SLR - 1 ) ) );
    fprintf( 1, ' z: %e\n', median( abs( m_z_CM / m_z_SLR - 1 ) ) );
    fprintf( 1, 'max( abs( CM - SLR ) )\n' );
    fprintf( 1, 'xy: %e\n', max( abs( m_xy_CM - m_xy_SLR ) ) );
    fprintf( 1, ' z: %e\n', max( abs( m_z_CM - m_z_SLR ) ) );

end
