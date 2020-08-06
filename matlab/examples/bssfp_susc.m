%% Susceptibility effect in balanced SSFP
% 
% This script numerically compares the configuration model with available
% analytical solutions for Lorentzian and spherical
% distributions.
% From a numerical point, the approach to determine the steady within the
% CM is very inefficient, but more accessible from a didactic point of view.

par = [];
opt = [];
str = [];

par.T1 = 150;
par.T2 = 120;
par.TR = 10;
par.flip_angle = 70;
par.T2p = 20;
par.sph_wdth = 0.9;
par.n_int = 25;
par.dummy = 10;

opt.T1 = [];
opt.T2 = [];
opt.TR = [];
opt.flip_angle = [];
opt.T2p = [];
opt.sph_wdth = [];
opt.n_int = [];
opt.dummy = [];

str.T1 = '[ms]';
str.T2 = '[ms]';
str.TR = '[ms]';
str.flip_angle = '[deg]';
str.T2p = 'T2'' (dephasing time for Lorentz distribution) [deg]';
str.sph_wdth = 'max. dephasing per TR [rad] (spherical distribution) divided by pi';
str.n_int = 'signal calculated at j * TR / n_int for 0 \le j \le n_int';
str.dummy = 'duration for equilibration [TR]';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    tau = par.TR / par.n_int;  % this is the highly ineffective part..
                   
    %% initialize configuration model
    
    % create instance
    
    cm = CoMoTk;
    
    % mandatory tissue parameters
    
    cm.R1 = 1 / par.T1;
    cm.R2 = 1 / par.T2;
    cm.D = 0;
        
    % RF parameters 
        
    RF_par = [];
    RF_par.FlipAngle = par.flip_angle * ( pi / 180.0 );
    RF_par.Phase = 0;                % non-alternating RF pulses
    
    % small time intervals
    
    Time_par = [];
    Time_par.lambda = 1;             % arbitrary unique index
    Time_par.tau = tau;

    %% calculate number of TR dummy cycles to establish a steady state
    
    n_TR = ceil( par.dummy * par.T1 / par.TR );

    % allocate space for results
    
    lorentz_cm = zeros( 1, par.n_int + 1 );
    spherical_cm = zeros( 1, par.n_int + 1 );
    x = linspace( 0, 1, par.n_int + 1 );
    tx = par.TR .* x;
        
    E2x_inv = exp( tx ./ par.T2 );
    
    for i = 1 : n_TR
        
        % execute RF pulse
        
        cm.RF( RF_par );                % RF pulse
        
        if ( i < n_TR )                 % dummy pulses 
               
            for j = 1 : par.n_int
                
                cm.time( Time_par );    % small time interval
                
            end
        
        else                            % here we get the steady-state result
            
            % parameters for readout
    
            sum_par = [];
            sum_par.omega = pi / par.TR;                     % since we did not use alternating RF pulses
            sum_par.R2p = 1 / par.T2p;                       % Lorentz distribution
            sum_par.bw_sphere = par.sph_wdth * pi / par.TR;  % spherical distribution
            
            for j = 0 : par.n_int
                
                if ( j > 0 )
                
                    cm.time( Time_par );    % small time interval
                
                end

                % Lorentz distribution
                
                cm.inhomogeneous_decay = @cm.lorentz_decay;
                res = cm.sum( sum_par );
                lorentz_cm( j + 1 ) = abs( res.xy );

                % spherical distribution
                
                cm.inhomogeneous_decay = @cm.spherical_decay;
                res = cm.sum( sum_par );
                spherical_cm( j + 1 ) = abs( res.xy );
               
                fprintf( 1, '%d / %d\n', j, par.n_int );

            end

        end    
        
        fprintf( 1, '%d / %d\n', i, n_TR );

    end
    
    % we ignore T2 decay to improve visibiliy of susceptibility effects 
    
    lorentz_cm = E2x_inv .* lorentz_cm;
    spherical_cm = E2x_inv .* spherical_cm;
    
    % calculate the analytical solution according to
    % Ganter, Magn Reson Med 2006, 56:687-691
    
    E1 = exp( - par.TR / par.T1 );
    E2 = exp( - par.TR / par.T2 );
    c = cosd( par.flip_angle );
    s = sind( par.flip_angle );
    a = 1 - E1 * c;
    b = c - E1;
    E22 = E2^2;
    Lambda = ( a - b * E22 - sqrt( ( a^2 - E22 * b^2 ) * ( 1 - E22 ) ) ) / ( a - b );
    A = 0.5 * ( a + b ) * E2 / ( a - 0.5 * Lambda * ( a - b ) );
    n = - n_TR : n_TR;
    tn = par.TR .* n;
    idx_0 = n_TR + 1;
    idx_neg = 1 : n_TR;
    m_n = A.^abs(n);
    m_n( idx_neg ) = - ( Lambda / ( E2 * A ) ) .* m_n( idx_neg );
    m_n = ( - 1i * ( 1 - E1 ) / ( a + b * Lambda ) * s ) .* m_n;
    m_n = ( -1 ).^n .* m_n;    % to simulate alternating RF pulses (or, alternatively, off-resonance)
    m_n = m_n( : );
    
    % calculate susceptibility dependent decay
    
    lorentz_decay = exp( - abs( tn' + tx ) ./ par.T2p );
    
    sigma = 0.5 * sum_par.bw_sphere;
    st_ = sigma .* ( tn' + tx );
    spherical_decay = besselj( 1, 2 .* st_ ) ./ ( st_ );
    spherical_decay( st_ == 0 ) = 1;
    
    lorentz_th = abs( sum( lorentz_decay .* m_n ) );
    spherical_th = abs( sum( spherical_decay .* m_n ) );
    
    plot( x, spherical_cm, 'b', x, spherical_th, 'b+', x, lorentz_cm, 'r', x, lorentz_th, 'r+' );
    
end
