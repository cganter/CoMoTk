%% Simple selective excitation
% This script shows how real RF pulses can be implemented.
% The example shows a simple amplitude modulated pulse
% The generalization to phase modulated and/or adiabatic pulses works along the same lines.
% Depending on the parameters, relaxation (e.g. for T2 < t_rf) and diffusion (large ADC) 
% effects may become relevant 

par = [];
opt = [];
str = [];

par.T1 = 100;
par.T2 = 10;
par.D = 0;
par.fa = 50;
par.sl_th = 1;
par.t_rf = 1;
par.supp_rf = 100;
par.qual_rf = 3;
par.filt_rf = 'Hamming';
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
    
    [ par, sel ] = sfv( par, opt, str );
    
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
    
    % rephasing gradient moment (for simplicity taken as - 0.5 * p_sl)
    % as discussed in the Handbook of MRI (and shown by the results),
    % this choice is reasonable but not optimal.
    
    p_refoc = [ 0; 0; - p_sl / 2 ];
    
    % assign unique handles for the two time intervals
    
    mu_rf = 1;
    mu_refoc = 2;
    
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
            param.mu = mu_rf;
            param.tau = tau_rf;
            param.p = p_rf;
            
            cm.time( param );
            
        else
            
            % in subsequent calls, we only need the handle
            
            param = [];
            param.mu = mu_rf;
            
            cm.time( param );
            
        end
        
        % and the rest of the small pulses
        
        param = [];
        param.FlipAngle = alpha( i + 1 );
        param.Phase = phase_rad;
        
        cm.RF( param );
        
    end
    
    %% collect the slice profile at locations x
    % prior to the rephasing gradient, the phase will vary across the slice
    
    m_iso_rf = zeros( 3, n_sl );
    
    for i = 1 : n_sl
        
        % we sum over all configurations (third argument == [])
        
        param = [];
        param.omega = omega;
        param.x = x( :, i );
        
        res = cm.sum( param );
        
        m_iso_rf( 1, i ) = real( res.xy );
        m_iso_rf( 2, i ) = imag( res.xy );
        m_iso_rf( 3, i ) = real( res.z );
        
    end
    
    %% apply the rephasing gradient
    % in absense of relaxation, the duration does not matter
    
    tau_refoc = 1;
    
    param = [];
    param.mu = mu_refoc;
    param.tau = tau_refoc;
    param.p = p_refoc;
    
    cm.time( param );
    
    %% collect the slice profile at locations x
    % now the phase should be essentially constant across the slice
    
    m_iso_refoc = zeros( 3, n_sl );
    
    for i = 1 : n_sl
        
        % a restriction of configurations is still not neccessary
        
        param = [];
        param.omega = omega;
        param.x = x( :, i );
        
        res = cm.sum( param );
        
        m_iso_refoc( 1, i ) = real( res.xy );
        m_iso_refoc( 2, i ) = imag( res.xy );
        m_iso_refoc( 3, i ) = real( res.z );
        
    end
    
    %% look at the results
    
    subplot( 1, 3, 1);
    plot( par.t_rf .* linspace( -0.5, 0.5, par.supp_rf + 1 ), alpha ./ max( alpha ) );
    xlim( [ - 0.5 * par.t_rf 0.5 * par.t_rf ] );
    title( 'B^+_1(t)' );
    
    subplot( 1, 3, 2 );
    plot( x( 3, : ), m_iso_rf( 1, : ), x( 3, : ), m_iso_rf( 2, : ), x( 3, : ), m_iso_rf( 3, : ) );
    legend( 'm_x', 'm_y', 'm_z' );
    title( 'after RF pulse' );
    
    subplot( 1, 3, 3 );
    plot( x( 3, : ), m_iso_refoc( 1, : ), x( 3, : ), m_iso_refoc( 2, : ), x( 3, : ), m_iso_refoc( 3, : ) );
    legend( 'm_x', 'm_y', 'm_z' );
    title( 'after rephasing gradient' );
    
end
