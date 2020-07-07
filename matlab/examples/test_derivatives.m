%% Check analytical derivatives
% The plots display the deviation between numerical and analytical derivatives:
% ( f( x + h ) - f( x ) ) / h - (df/dx)( x )
% The result is plotted logarithmically against h over three decades.
% If the analytical derivative is calculated correctly, we expect
% the deviation to be proportional to h (in leading order)

% test is based on individual configurations or isochromats

par = [];
opt = [];
str = [];

par.derivative = 'R1';
par.check = 'isochromats';
par.tissues = 3;
par.diffusion = 'none';
par.coupling = 'none';
par.shape = 'implicit';
par.train = 50;
par.d = 2;
par.n_fa = 3;
par.n_ph = 3;

opt.derivative = { 'R1', 'R2', 'D', 'B1', 'FlipAngle', 'Phase', 'tau', 'p', 's', 'S' };
opt.check = { 'configurations', 'isochromats' };
opt.tissues = [];
opt.diffusion = { 'none', 'isotropic', 'tensor' };
opt.coupling = { 'none', 'magnetization transfer', 'J-coupling' };
opt.shape = { 'implicit', 'explicit' };
opt.train = [];
opt.d = [];
opt.n_fa = [];
opt.n_ph = [];

str.derivative = 'variable parameter';
str.check = 'which output to check';
str.tissues = 'number of tissues';
str.diffusion = 'how to include diffusion effects';
str.coupling = 'magnetization transfer/exchange, J-coupling';
str.shape = 'gradient shape';
str.train = 'number of intervals, after which the derivatives are checked';
str.d = 'dimension of configuration model == number of different intervals (< train)';
str.n_fa = 'number of different flip angles (< train)';
str.n_ph = 'number of different phases (< train)';

% typical scale of parameters

scale = struct( ...
    'R1', 1e-2, ...
    'R2', 1e-1, ...
    'D', 1, ...
    'B1', 1, ...
    'FlipAngle', pi, ...
    'Phase', 2 * pi, ...
    'tau', 1, ...
    'p', 1e-2, ...
    's', 1, ...
    'S', 1e-2 ...
    );

% parameter degrees of freedom

deg = struct( ...
    'R1', 1, ...
    'R2', 1, ...
    'D', 1, ...
    'B1', 1, ...
    'FlipAngle', 1, ...
    'Phase', 1, ...
    'tau', 1, ...
    'p', 3, ...
    's', 3, ...
    'S', 1 ...
    );

% logarithmic spacing over three decades (1e-5 ... 1e-2)

x = 1e-5 * ( 1e3.^( 0 : 0.1 : 1 ) );
n_x = length( x );
    
while ( true )
    
    [ par, sel ] = sfv( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end

    % rule out invalid combinations
    
    if ( isequal( par.diffusion, 'none' ) && ...
        ( isequal( par.derivative, 'p' ) || isequal( par.derivative, 's' ) || isequal( par.derivative, 'S' ) || isequal( par.shape, 'explicit' ) ) )
        
        fprintf( 1, 'Gradient moment irrelevant without diffusion.\n\n' );
        continue;
        
    end
    
    if ( isequal( par.diffusion, 'none' ) && isequal( par.derivative, 'D' ) )
        
        fprintf( 1, 'Diffusion derivative requires diffusion...\n\n' );
        continue;
        
    end
    
    if ( ( isequal( par.derivative, 'p' ) || isequal( par.derivative, 's' ) || isequal( par.derivative, 'S' ) ) && isequal( par.shape, 'implicit' ) )
    
        fprintf( 1, 'Derivative with respect to shape requires the latter to be explicit.\n' );
        fprintf( 1, 'This also includes derivatives with respect to p.\n\n' );
        continue;
        
    end

    if ( ~isequal( par.coupling, 'none' ) )
    
        fprintf( 1, 'Spin coupling not implemented yet.\n' );
        continue;
        
    end

    % set all sequence parameters
    
    var_.R1 = scale.R1 .* ( 0.5 + rand( 1, par.tissues ) );
    var_.R2 = scale.R2 .* ( 0.5 + rand( 1, par.tissues ) );
    
    if ( isequal( par.diffusion, 'isotropic' ) )
        
        var_.D = scale.D .* ( 0.5 + rand( 1, par.tissues ) );
        deg.D = 1;
        
    elseif ( isequal( par.diffusion, 'tensor' ) )
        
        var_.D = zeros( 3, 3, par.tissues );
        
        for i = 1 : par.tissues
            
            sqrt_D = randn( 3 );
            var_.D( :, :, i ) = sqrt_D * sqrt_D';
            var_.D( :, :, i ) = ( 3 * scale.D / trace( var_.D( :, :, i ) ) ) .* var_.D( :, :, i );

            deg.D = 6;
            
        end
        
    elseif ( isequal( par.diffusion, 'none' ) )
        
        var_.D = zeros( 1, par.tissues );

        deg.D = 1;
        
    end

    var_.B1 = 1;

    var_.FlipAngle = scale.FlipAngle .* rand( 1, par.n_fa );
    var_.Phase = scale.Phase .* rand( 1, par.n_ph );
    var_.tau = scale.tau .* ( 0.5 + rand( 1, par.d ) );
    var_.p = scale.p .* randn( 3, par.d );
    
    if ( isequal( par.shape, 'explicit' ) )
        
        var_.s = scale.s .* randn( 3, par.d );
        
        if ( isequal( par.diffusion, 'isotropic' ) )
            
            var_.S = scale.S .* ( 0.5 + rand( 1, par.d ) );
            
            deg.S = 1;
            
        elseif ( isequal( par.diffusion, 'tensor' ) )
            
            var_.S = zeros( 3, 3, par.d );
            
            for i = 1 : par.d
                
                sqrt_S = randn( 3 );
                var_.S( :, :, i ) = sqrt_S * sqrt_S';
                var_.S( :, :, i ) = ( 3 * scale.S / trace( var_.S( :, :, i ) ) ) .* var_.S( :, :, i );

                deg.S = 6;
                
            end

        end
            
    end

    % define the variable increment
    
    if ( ismember( par.derivative, { 'R1', 'R2', 'B1', 'FlipAngle', 'Phase', 'tau' } ) )
        
        dx = scale.( par.derivative );
        dx_vec = dx;
        
    elseif ( ismember( par.derivative, { 'p', 's' } ) )
        
        dx = scale.( par.derivative ) .* randn( 3, 1 );
        dx_vec = dx;
        
    elseif ( ismember( par.derivative, { 'D', 'S' } ) )
        
        if ( isequal( par.diffusion, 'isotropic' ) )
            
            dx = scale.( par.derivative );
            dx_vec = dx;
            
        elseif ( isequal( par.diffusion, 'tensor' ) )

            sqrt_dx = randn( 3 );
            dx = sqrt_dx * sqrt_dx';
            dx = ( 3 * scale.( par.derivative ) / trace( dx ) ) .* dx;
            dx_vec = [ dx( 1, 1 ); dx( 2, 2 ); dx( 3, 3 ); dx( 1, 2 ); dx( 1, 3 ); dx( 2, 3 ) ];
            
        end
        
    end
            
    % randomized sequence of RF pulses and time intervals
    
    idx.fa = randi( par.n_fa, 1, par.train );
    idx.ph = randi( par.n_ph, 1, par.train );
    idx.tau = randi( par.d, 1, par.train );
    
    % randomly select derivatives to calculate
        
    if ( par.tissues > 1 && ismember( par.derivative, { 'R1', 'R2', 'D' } ) )
        
        idx.der = randperm( par.tissues );
        idx.der = idx.der( 1 : 2 );
        
    elseif ( par.d > 1 && ismember( par.derivative, { 'tau', 'p', 's', 'S' } ) )
        
        idx.der = randperm( par.d );
        idx.der = idx.der( 1 : 2 );
                
    elseif ( par.n_fa > 1 && isequal( par.derivative, { 'FlipAngle' } ) )
        
        idx.der = randperm( par.n_fa );
        idx.der = idx.der( 1 : 2 );
                
    elseif ( par.n_ph > 1 && isequal( par.derivative, { 'Phase' } ) )
        
        idx.der = randperm( par.n_ph );
        idx.der = idx.der( 1 : 2 );
                
    else
        
        idx.der = 1;
        
    end
    
    disp( idx.der );
    
    % balanced sequences are most sensitive on the pure shape S

    if ( isequal( par.derivative, 'S' ) )
        
        var_.p( : ) = 0;
        
    end
        
    % do tests
    
    fprintf( 1, 'Check derivative with respect to %s.\n', par.derivative );
        
    rel_diff = zeros( 1, n_x );
    var = var_;
    
    % calculate the for various arguments
    
    for i_x = 0 : n_x
        
        fprintf( 1, '%d / %d\n', i_x, n_x );
        
        % set variable parameter
   
        if ( i_x > 0 )
            
            if ( ismember( par.derivative, { 'R1', 'R2', 'B1', 'FlipAngle', 'Phase', 'tau' } ) )
                
                var.( par.derivative )( idx.der( 1 ) ) = var_.( par.derivative )( idx.der( 1 ) ) + x( i_x ) * dx;

            elseif ( ismember( par.derivative, { 'p', 's' } ) )
                
                var.( par.derivative )( :, idx.der( 1 ) ) = var_.( par.derivative )( :, idx.der( 1 ) ) + x( i_x ) * dx;
                
            elseif ( ismember( par.derivative, { 'D', 'S' } ) )
                
                if ( isequal( par.diffusion, 'isotropic' ) )
                    
                    var.( par.derivative )( idx.der( 1 ) ) = var_.( par.derivative )( idx.der( 1 ) ) + x( i_x ) * dx;

                elseif ( isequal( par.diffusion, 'tensor' ) )
                    
                    var.( par.derivative )( :, :, idx.der( 1 ) ) = var_.( par.derivative )( :, :, idx.der( 1 ) ) + x( i_x ) * dx;

                end
                
            end
            
        end
        
        % initialize configuration model
                        
        cm = CoMoTk;

        cm.R1 = var.R1;
        cm.R2 = var.R2;
        cm.D = var.D;
        cm.B1 = var.B1;
        
        % get default options
        
        options = cm.options;
        
        options.alloc_n = 5000;
        options.alloc_d = par.d;
        options.epsilon = 0;
        
        % set new options
        
        cm.options = options;
        
        % initialize everything
        
        cm.init_configuration ( [ zeros( 1, par.tissues ); zeros( 1, par.tissues ); ones( 1, par.tissues ) ] );
        
        % specify derivatives to calculate
        
        param = [];
        param.( par.derivative ) = idx.der;

        cm.set_derivatives( param );

        % execute the sequence
        
        for i_n = 1 : par.train
            
            % RF pulse
            
            param = [];
            param.FlipAngle = var.FlipAngle( idx.fa( i_n ) );
            param.Phase = var.Phase( idx.ph( i_n ) );
            param.handle_FlipAngle = idx.fa( i_n );
            param.handle_Phase = idx.ph( i_n );
            
            cm.RF( param );

            % time interval
            
            param = [];
            param.mu = idx.tau( i_n );
            param.tau = var.tau( idx.tau( i_n ) );
            param.p = var.p( :, idx.tau( i_n ) );
            
            if ( isequal( par.shape, 'explicit' ) )
                
                param.s = var.s( :, idx.tau( i_n ) );
                
                if ( isequal( par.diffusion, 'isotropic' ) )
                    
                    param.S = var.S( idx.tau( i_n ) );
                    
                elseif ( isequal( par.diffusion, 'tensor' ) )
                
                    param.S = var.S( :, :, idx.tau( i_n ) );

                end
                
            end
            
            cm.time( param );
            
        end
        
        if ( isequal( par.check, 'isochromats' ) )
            
            if ( i_x == 0 )
                
                param = [];
                
                res0 = cm.sum( param );
                diso0.dxy = res0.dxy.( [ 'd', par.derivative ] )( :, 1 );
                diso0.dz = res0.dz.( [ 'd', par.derivative ] )( :, 1 );
                lin0 = sum( [ diso0.dxy, diso0.dz ] .* dx_vec, 1 );
                
                disp( lin0 );
                
            else
                
                param = [];
                
                res = cm.sum( param );
                
                d_m = ( [ res.xy - res0.xy, res.z - res0.z ] ) ./ x( i_x );
                       
                disp( d_m );
                
                d_m = d_m - lin0;
                
            end
            
        elseif ( isequal( par.check, 'configurations' ) )

            if ( i_x == 0 )
                
                m0 = cm.m( :, :, cm.b_n );
                d_m0 = cm.dm.( [ 'd', par.derivative ] )( :, :, cm.b_n, :, 1 );
                
            else
                
                d_m = ( cm.m( :, :, cm.b_n ) - m0 ) ./ x( i_x ) - ...
                    sum( d_m0 .* reshape( dx_vec, [ 1, 1, 1, length( dx_vec ) ] ), 4 );
                
            end
                        
        end
        
        if ( i_x > 0 )
            
            rel_diff( i_x ) = sqrt( real( d_m( : )' * d_m( : ) ) );
            
        end
        
    end
    
    subplot( 1, 1, 1 );
    
    loglog( x, rel_diff, '*', x, rel_diff( 1 ) .* x ./ x( 1 ) );
    title( par.derivative );
    
end
