%% Check analytical derivatives
% The plots display the deviation between numerical and analytical derivatives:
% ( f( x + h ) - f( x ) ) / h - (df/dx)( x )
% The result is plotted logarithmically against h over three decades.
% If the analytical derivative is calculated correctly, we expect
% the deviation to be proportional to h (in leading order)

% test is based on individual configurations or isochromats

b_iso = [];

while( isempty( b_iso ) )
    
    a_ = input( 'Check with configurations (1) or isochromats (2)?\n' );

    if ( a_ == 1 )
        
        b_iso = false;
        
        fprintf( 1, 'Check with configurations\n' );
        
    elseif ( a_ == 2 )
        
        b_iso = true;
        
        fprintf( 1, 'Check with isochromats\n' );
        
    else
        
        fprintf( '1 or 2 expected.\' );
        
    end
    
end

% test should be run +/- diffusion effects

b_diff = [];

while( isempty( b_diff ) )
    
    a_ = input( 'Calculate with diffusion? (1 == yes, 2 == no)\n' );

    if ( a_ == 1 )
        
        b_diff = true;
        
        fprintf( 1, 'Check derivatives with diffusion effects\n' );
        
    elseif ( a_ == 2 )
        
        b_diff = false;
        
        fprintf( 1, 'Check derivatives without diffusion effects\n' );
        
    else
        
        fprintf( '1 or 2 expected.\' );
        
    end
    
end

% in case of diffusion, the shape can be specified explicitly or not
% (the latter case assumes a constant gradient in every time period)

if ( b_diff )
    
    b_shape = [];
    
    while( isempty( b_shape ) )
        
        a_ = input( 'Calculate with explicit shape? (1 == yes, 2 == no)\n' );
        
        if ( a_ == 1 )
            
            b_shape = true;
            
            fprintf( 1, 'Variable gradient shape\n' );
            
        elseif ( a_ == 2 )
            
            b_shape = false;
            
            fprintf( 1, 'Constant gradient\n' );
            
        else
            
            fprintf( '1 or 2 expected.\' );
            
        end
        
    end
        
end

% dimensions

d = 2;    % instantaneous RF pulses and two different time intervals
n_t = 3;  % number of tissues
n_al = 3; % number of different flip angles
n_ph = 3; % number of different phases

% typical scale of parameters

R1_scale = 1e-3;
R2_scale = 1e-2;
D_scale = 1;
B1_scale = 1;
FlipAngle_scale = pi;
Phase_scale = 2 * pi;
tau_scale = 1;
p_scale = [ 1e-2; 1e-2; 1e-2 ];
s_scale = [ 1e-2; 1e-2; 1e-2; 1e-4 ];

% logarithmic spacing over three decades

x = 1e-5 * ( 1000.^( 0 : 0.1 : 1 ) );

n_x = length( x );

% number of intervals, after which we check the derivatives

n = 50;

% random settings

idx_alpha = ceil( n_al .* rand( 1, n ) );
idx_phase = ceil( n_ph .* rand( 1, n ) );
idx_tau = ceil( d .* rand( 1, n ) );

% randomly define a subset of derivatives to calculate

idx_ = randperm( n_t );
idx_dR1 = idx_( 1 : 2 );

idx_ = randperm( n_t );
idx_dR2 = idx_( 1 : 2 );

idx_ = randperm( n_t );
idx_dD = idx_( 1 : 2 );

idx_ = randperm( n_al );
idx_dFlipAngle = idx_( 1 : 2 );

idx_ = randperm( n_ph );
idx_dPhase = idx_( 1 : 2 );

idx_ = randperm( d );
idx_dtau = idx_( 1 : 2 );

idx_ = randperm( d );
idx_dp = idx_( 1 : 2 );

idx_ = randperm( d );
idx_ds = idx_( 1 : 2 );

% derivative strings

if ( b_diff )

    if ( b_shape )
        
        der_str = { 'R1', 'R2', 'D', 'B1', 'FlipAngle', 'Phase', 'tau', 'p', 's' };

    else
        
        der_str = { 'R1', 'R2', 'D', 'B1', 'FlipAngle', 'Phase', 'tau', 'p' };
        
    end
        
else

    der_str = { 'R1', 'R2', 'B1', 'FlipAngle', 'Phase', 'tau' };

end

% do tests

for i_d = 1 : length( der_str )
    
    fprintf( 1, 'Check derivative with respect to %s.\n', der_str{ i_d } );
       
    % set various parameters
    
    R1 = R1_scale .* ( 1 + rand( 1, n_t ) );
    R2 = R2_scale .* ( 1 + rand( 1, n_t ) );
    
    if ( b_diff )
        
        D = D_scale .* ( 1 + rand( 1, n_t ) );
        
    else
        
        D = zeros( 1, n_t );
        
    end
  
    B1 = 1;
    alpha_0 = FlipAngle_scale .* rand( 1, n_al );
    phase_0 = Phase_scale .* rand( 1, n_ph );
    tau_0 = tau_scale .* ( 1 + rand( 1, d ) );
    p_0 = p_scale .* randn( 3, d );
    s_0 = s_scale .* randn( 4, d );
    
    % random gradient orientations
    
    e3 = randn( 3, 1 );
    e3 = e3 ./ sqrt( e3' * e3 );
    dp = p_scale .* e3;
    
    e4 = randn( 4, 1 );
    e4 = e4 ./ sqrt( e4' * e4 );
    ds = s_scale .* e4;
    
    rel_diff = zeros( 1, n_x );

    if ( isequal( der_str{ i_d }, 'R1' ) )
        
        idx_d = idx_dR1;
        
    elseif ( isequal( der_str{ i_d }, 'R2' ) )
        
        idx_d = idx_dR2;
        
    elseif ( isequal( der_str{ i_d }, 'D' ) )
        
        idx_d = idx_dD;
        
    elseif ( isequal( der_str{ i_d }, 'B1' ) )
        
        idx_d = 1;
        
    elseif ( isequal( der_str{ i_d }, 'FlipAngle' ) )
        
        idx_d = idx_dFlipAngle;
        
    elseif ( isequal( der_str{ i_d }, 'Phase' ) )
        
        idx_d = idx_dPhase;
        
    elseif ( isequal( der_str{ i_d }, 'tau' ) )
        
        idx_d = idx_dtau;
        
    elseif ( isequal( der_str{ i_d }, 'p' ) )
        
        idx_d = idx_dp;
        
    elseif ( isequal( der_str{ i_d }, 's' ) )
        
        idx_d = idx_ds;
        
    else
        
        error( 'Unknown parameter.' );
        
    end

    for i_x = 0 : n_x
        
        fprintf( 1, '%d / %d\n', i_x, n_x );
        
        % initialize configuration model
        
        cm = CoMoTk;
        
        % parameters
        
        cm.R1 = R1;
        cm.R2 = R2;
        cm.D = D;
        cm.B1 = B1;

        alpha = alpha_0;
        phase = phase_0;
        tau = tau_0;
        p = p_0;
        s = s_0;
        
        % increment on the variable
        
        if ( i_x > 0 )
            
            if ( isequal( der_str{ i_d }, 'R1' ) )
                
                cm.R1( idx_dR1( 1 ) ) = R1( idx_dR1( 1 ) ) + x( i_x ) * R1_scale;
                
            elseif ( isequal( der_str{ i_d }, 'R2' ) )
                
                cm.R2( idx_dR2( 1 ) ) = R2( idx_dR2( 1 ) ) + x( i_x ) * R2_scale;
                
            elseif ( isequal( der_str{ i_d }, 'D' ) )
                
                cm.D( idx_dD( 1 ) ) = D( idx_dD( 1 ) ) + x( i_x ) * D_scale;
                
            elseif ( isequal( der_str{ i_d }, 'B1' ) )
                
                cm.B1 = B1 + x( i_x ) * B1_scale;
                
            elseif ( isequal( der_str{ i_d }, 'FlipAngle' ) )
                
                alpha( idx_dFlipAngle( 1 ) ) = alpha_0( idx_dFlipAngle( 1 ) ) + x( i_x ) * FlipAngle_scale;
                
            elseif ( isequal( der_str{ i_d }, 'Phase' ) )
                
                phase( idx_dPhase( 1 ) ) = phase_0( idx_dPhase( 1 ) ) + x( i_x ) * Phase_scale;
                
            elseif ( isequal( der_str{ i_d }, 'tau' ) )
                
                tau( idx_dtau( 1 ) ) = tau_0( idx_dtau( 1 ) ) + x( i_x ) * tau_scale;
                
            elseif ( isequal( der_str{ i_d }, 'p' ) )
                
                p( :, idx_dp( 1 ) ) = p_0( :, idx_dp( 1 ) ) + x( i_x ) * dp;
                
            elseif ( isequal( der_str{ i_d }, 's' ) )
                
                s( :, idx_ds( 1 ) ) = s_0( :, idx_ds( 1 ) ) + x( i_x ) * ds;
                               
            end
            
        end
        
        % get default options
        
        options = cm.options;
        
        options.alloc_n = 5000;
        options.alloc_d = d;
        options.epsilon = 0;
        
        % set new options
        
        cm.options = options;
        
        % initialize everything
        
        cm.init_configuration ( [ zeros( 1, n_t ); zeros( 1, n_t ); ones( 1, n_t ) ] );
        
        % specify derivatives to calulate

        if ( ~isequal( der_str{ i_d }, 'B1' ) )
        
            cm.set_derivatives ( der_str{ i_d }, idx_d );
            
        else
            
            cm.set_derivatives ( der_str{ i_d } );

        end
                    
        for i_n = 1 : n
                        
            % RF pulse
                        
            cm.RF( alpha( idx_alpha( i_n ) ), phase( idx_phase( i_n ) ) , 'FlipAngle', idx_alpha( i_n ) , 'Phase', idx_phase( i_n ) );
            
            % time interval

            if ( b_diff && b_shape )
            
                cm.time( idx_tau( i_n ), 'tau', tau( idx_tau( i_n ) ), 'p', p( :, idx_tau( i_n ) ) , 's', s( :, idx_tau( i_n ) ) );
                
            else
                
                cm.time( idx_tau( i_n ), 'tau', tau( idx_tau( i_n ) ), 'p', p( :, idx_tau( i_n ) ) );
                
            end
                       
        end

        if ( b_iso )
            
            if ( i_x == 0 )
                
                iso0 = cm.isochromat( 0, [], [] );
                
            else
                
                iso = cm.isochromat( 0, [], [] );
                
            end
            
        else
            
            if ( isequal( der_str{ i_d }, 'R1' ) )
                
                if ( i_x == 0 )
                    
                    
                    m0 = cm.m( :, idx_dR1( 1 ), cm.b_n );
                    d_m0 = cm.dm_dR1( :, 1, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, idx_dR1( 1 ), cm.b_n ) - m0 ) ./ ( cm.R1( idx_dR1( 1 ) ) - R1( idx_dR1( 1 ) ) ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'R2' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, idx_dR2( 1 ), cm.b_n );
                    d_m0 = cm.dm_dR2( :, 1, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, idx_dR2( 1 ), cm.b_n ) - m0 ) ./ ( cm.R2( idx_dR2( 1 ) ) - R2( idx_dR2( 1 ) ) ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'D' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, idx_dD( 1 ), cm.b_n );
                    d_m0 = cm.dm_dD( :, 1, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, idx_dD( 1 ), cm.b_n ) - m0 ) ./ ( cm.D( idx_dD( 1 ) ) - D( idx_dD( 1 ) ) ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'B1' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, :, cm.b_n );
                    d_m0 = cm.dm_dB1( :, :, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ ( cm.B1 - B1 ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'FlipAngle' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, :, cm.b_n );
                    d_m0 = cm.dm_dFlipAngle( :, :, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ ( alpha( idx_dFlipAngle( 1 ) ) - alpha_0( idx_dFlipAngle( 1 ) ) ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'Phase' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, :, cm.b_n );
                    d_m0 = cm.dm_dPhase( :, :, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ ( phase( idx_dPhase( 1 ) ) - phase_0( idx_dPhase( 1 ) ) ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'tau' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, :, cm.b_n );
                    d_m0 = cm.dm_dtau( :, :, cm.b_n );
                    
                else
                    
                    d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ ( tau( idx_dtau( 1 ) ) - tau_0( idx_dtau( 1 ) ) ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 'p' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, :, cm.b_n );
                    d_m0 = sum( cm.dm_dp( :, :, cm.b_n, :, 1 ) .* reshape( dp, [ 1, 1, 1, 3 ] ), 4 );
                    
                else
                    
                    d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ x( i_x ) - d_m0;
                    
                end
                
            elseif ( isequal( der_str{ i_d }, 's' ) )
                
                if ( i_x == 0 )
                    
                    m0 = cm.m( :, :, cm.b_n );
                    d_m0 = sum( cm.dm_ds( :, :, cm.b_n, :, 1 ) .* reshape( ds, [ 1, 1, 1, 4 ] ), 4 );
                    
                else
                    
                    d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ x( i_x ) - d_m0;
                    
                end
                
            end
            
        end
        
        if ( i_x > 0 )
            
            if( b_iso )
                
                if ( isequal( der_str{ i_d }, 'R1' ) )
                    
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( cm.R1( idx_dR1( 1 ) ) - R1( idx_dR1( 1 ) ) ) - [ iso0.dm_dR1.xy( 1 ), iso0.dm_dR1.z( 1 ) ];
                    
                elseif ( isequal( der_str{ i_d }, 'R2' ) )
                
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( cm.R2( idx_dR2( 1 ) ) - R2( idx_dR2( 1 ) ) ) - [ iso0.dm_dR2.xy( 1 ), iso0.dm_dR2.z( 1 ) ];
                    
                elseif ( isequal( der_str{ i_d }, 'D' ) )
                
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( cm.D( idx_dD( 1 ) ) - D( idx_dD( 1 ) ) ) - [ iso0.dm_dD.xy( 1 ), iso0.dm_dD.z( 1 ) ];
                    
                elseif ( isequal( der_str{ i_d }, 'B1' ) )
                
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( cm.B1 - B1 ) - [ iso0.dm_dB1.xy, iso0.dm_dB1.z ];
                    
                elseif ( isequal( der_str{ i_d }, 'FlipAngle' ) )
                
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( alpha( idx_dFlipAngle( 1 ) ) - alpha_0( idx_dFlipAngle( 1 ) ) ) - [ iso0.dm_dFlipAngle.xy( 1 ), iso0.dm_dFlipAngle.z( 1 ) ];
                    
                elseif ( isequal( der_str{ i_d }, 'Phase' ) )
                
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( phase( idx_dPhase( 1 ) ) - phase_0( idx_dPhase( 1 ) ) ) - [ iso0.dm_dPhase.xy( 1 ), iso0.dm_dPhase.z( 1 ) ];
                    
                elseif ( isequal( der_str{ i_d }, 'tau' ) )
                
                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ ( tau( idx_dtau( 1 ) ) - tau_0( idx_dtau( 1 ) ) ) - [ iso0.dm_dtau.xy( 1 ), iso0.dm_dtau.z( 1 ) ];
                    
               elseif ( isequal( der_str{ i_d }, 'p' ) )

                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ x( i_x ) - ...
                        sum( [ iso0.dm_dp.xy( :, 1 ), iso0.dm_dp.z( :, 1 ) ] .* dp );
                    
               elseif ( isequal( der_str{ i_d }, 's' ) )

                    d_m = ( [ iso.xy - iso0.xy, iso.z - iso0.z ] ) ./ x( i_x ) - ...
                        sum( [ iso0.dm_ds.xy( :, 1 ), iso0.dm_ds.z( :, 1 ) ] .* ds );
                    
                end
                    
            end
            
            rel_diff( i_x ) = sqrt( real( d_m( : )' * d_m( : ) ) );
            
        end
        
    end
        
    subplot( ceil( length( der_str ) / 3 ), 3, i_d );
    
    loglog( x, rel_diff, '*', x, rel_diff( end ) .* x ./ x( end ) );
    title( der_str{ i_d } );
 
end
