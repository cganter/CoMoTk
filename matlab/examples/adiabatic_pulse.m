%% Simple selective excitation
% This script shows how real RF pulses can be implemented.
% The example shows a simple amplitude modulated pulse
% The generalization to phase modulated and/or adiabatic pulses works along the same lines.
% Depending on the parameters, relaxation (e.g. for T2 < t_rf) and diffusion (large ADC) 
% effects may become relevant 

% fixed parameters from Garwood et al., JMR (2001) 153, 155-177
% (Table 1 therein)

A = 2 * pi * 25; % [rad kHz] half RF pulse bandwidth 
t_rf = 2;        % [ms] RF pulse duration

% HS1 pulse

beta = asech( 0.01 );
mu = A / beta;
om1_HS1 = 2 * pi * 7.56;   % [rad kHz] rotation frequency associated with peak B1+

% Chirp pulse

om1_Chirp = 2 * pi * 3.38;   % [rad kHz] rotation frequency associated with peak B1+

% displayed number of freqencies

n_freq = 201;

% free parameters

par = [];
opt = [];
str = [];

par.T1 = Inf;
par.T2 = 1;
par.D = 0;
par.n_rf = 500;

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.n_rf = [];

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.n_rf = 'number of RF support points';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    % time interval between RF pulses
    
    n_rf = par.n_rf; 
    dt = t_rf / ( n_rf - 1 );
    t = linspace( -1, 1, n_rf );

    % off-resonance frequency [rad kHz]
    
    Omega = linspace( -A, A, n_freq );
    
    % define adiabatic pulse
    
    alpha_HS1 = ( om1_HS1 * dt ) .* sech( beta .* t );
    phi_HS1 = mu .* log( sech( beta .* t ) );
        
    alpha_Chirp = ( om1_Chirp * dt ) .* ones( size( t ) );
    phi_Chirp = - ( 0.5 * A ) .* t.^2;
        
    %% initialize configuration model
    
    cm_HS1 = CoMo;
    cm_HS1_0 = CoMo;
    cm_Chirp_0 = CoMo;
    
    % mandatory tissue parameters
    
    cm_HS1.R1 = 1 / par.T1;
    cm_HS1.R2 = 1 / par.T2;
    cm_HS1.D = par.D;
    
    cm_HS1_0.R1 = 0;
    cm_HS1_0.R2 = 0;
    cm_HS1_0.D = 0;
    
    cm_Chirp_0.R1 = 0;
    cm_Chirp_0.R2 = 0;
    cm_Chirp_0.D = 0;
    
    % configuration space resolution
    
    cm_HS1.d_tau = dt;
    cm_HS1_0.d_tau = dt;
    cm_Chirp_0.d_tau = dt;

    % allocate space for transverse and longitudinal magnetization during the RF pulse
   
    mz_HS1 = zeros( 1, n_rf );
    mz_HS1_rapid = zeros( 1, n_rf );
    mz_Chirp_rapid = zeros( 1, n_rf );
    
    mz_freq_HS1 = zeros( 1, n_freq );
    mz_freq_HS1_rapid = zeros( 1, n_freq );
    mz_freq_Chirp_rapid = zeros( 1, n_freq );
    
    %% execute the RF pulse
    
    for i = 1 : n_rf

        % small RF pulse
            
        param = [];
        param.FlipAngle = alpha_HS1( i );
        param.Phase = phi_HS1( i );

        cm_HS1.RF( param );  
        cm_HS1_0.RF( param );  
        
        param = [];
        param.FlipAngle = alpha_Chirp( i );
        param.Phase = phi_Chirp( i );
        
        cm_Chirp_0.RF( param );            
            
        % get actual state
        % (== complete sum over configurations)
        
        param = [];

        res = cm_HS1.sum( param );
        mz_HS1( i ) = real( res.z );

        res = cm_HS1_0.sum( param );
        mz_HS1_rapid( i ) = real( res.z );

        res = cm_Chirp_0.sum( param );
        mz_Chirp_rapid( i ) = real( res.z );

        if ( i == n_rf )
            
            for j = 1 : n_freq
                
                param.omega = Omega( j );        
                
                res = cm_HS1.sum( param );
                mz_freq_HS1( j ) = real( res.z );
                
                res = cm_HS1_0.sum( param );
                mz_freq_HS1_rapid( j ) = real( res.z );
                
                res = cm_Chirp_0.sum( param );
                mz_freq_Chirp_rapid( j ) = real( res.z );
                
            end
        
        end
        
        % time interval
            
        param = [];
        param.tau = dt;
           
        cm_HS1.time( param );
        cm_HS1_0.time( param );
        cm_Chirp_0.time( param );
            
    end
    
    %% look at the results

    t = ( 0.5 * t_rf ) .* ( t + 1 );
    
    ax = subplot( 1, 2, 1 );
    plot( t, mz_HS1_rapid, '-b', t, mz_Chirp_rapid, '-r', t, mz_HS1, '--b' );
    
    xlabel( '$t$ [ms]', 'Interpreter', 'latex' );
    ylabel( '$m_z\left(t\right)$', 'Interpreter', 'latex' );
    title( 'Adiabatic Transient Phase', 'Interpreter', 'latex' );
    
    width = 14;
    height = 7.5;

    set( gcf, 'Units', 'centimeters' );
    set( gcf, 'Position', [ 0, 0, width, height ] );
    set( gcf, 'Color', 'w' );

    delta = 0.01;
    
    ti = ax.TightInset;
    left = ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.5 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    ax = subplot( 1, 2, 2 );
    plot( Omega ./ ( 2 * pi ), mz_freq_HS1_rapid, '-b', Omega ./ ( 2 * pi ), mz_freq_Chirp_rapid, '-r', Omega ./ ( 2 * pi ), mz_freq_HS1, '--b' );
        
    xlabel( '$\omega/2\pi$ [kHz]', 'Interpreter', 'latex' );
    ylabel( '$m_z\left(\omega\right)$', 'Interpreter', 'latex' );
    title( 'Inversion Profile', 'Interpreter', 'latex' );
    legend( 'HS1, $T_2 = \infty$', 'chirp, $T_2 = \infty$', [ 'HS1, $T_2 = $ ', num2str( par.T2 ) ], 'Interpreter', 'latex', 'Location', 'north' );

    ti = ax.TightInset;
    left = 0.5 + ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 0.5 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];

end
