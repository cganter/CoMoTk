%% Compares diffusion of unbalanced SSFP with available analytical solution 
% derived by D. Freed et al. in: J Chem Phys (2001) 115, 4249-4258

%% Simulation parameters
% (by default, all set by the user)

par = [];
opt = [];
str = [];

par.T1 = 500;
par.T2 = 50;
par.D = 3;
par.TR = 5;
par.fa = 30;
par.dx_2pi = 20;
par.n_max = 8;
par.t_prep = 10;
par.tensor = 'no';
par.implicit = 'no';

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.TR = [];
opt.fa = [];
opt.dx_2pi = [];
opt.n_max = [];
opt.t_prep = [];
opt.tensor = { 'yes', 'no' };
opt.implicit = { 'yes', 'no' };

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/ms] isotropic ADC';
str.TR = '[ms] repetition time';
str.fa = '[deg] flip angle';
str.dx_2pi = '[um] scale on which crusher effects dephasing of 2 * pi'; 
str.n_max = 'maximum configuration order';
str.t_prep = 'duration of preparation phase in units of T1';
str.tensor = 'test diffusion tensor implementation';
str.implicit = 'test implicit vs. explicit gradient shape implementation';

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    % read options
    
    b_tensor = false;
    b_implicit = false;
    
    if ( isequal( par.tensor, 'yes' ) )
        
        b_tensor = true;
        
    end
        
    if ( isequal( par.implicit, 'yes' ) )
        
        b_implicit = true;
        
    end
        
    % number of TR cycles to approach steady state
    
    num_TR = ceil( par.t_prep * par.T1 / par.TR );
    
    % convert to units, as expected by CoMoTk
    
    fa_rad = par.fa * pi / 180;
    
    % routine returns configuration orders n = - n_max : n_max - 1
    
    n_max = par.n_max;
    
    % random orthogonal matrix
    % (== rotation by random angle 'al_' around random axis 'no_')
    
    al_ = 2 * pi * rand( 1 );

    c_al = cos( al_ );
    s_al = sin( al_ );
    
    ph_ = 2 * pi * rand( 1 );

    c_ph = cos( ph_ );
    s_ph = sin( ph_ );

    th_ = pi * rand( 1 );
    
    c_th = cos( th_ );
    s_th = sin( th_ );
        
    no_ = [ s_th * c_ph; s_th * s_ph; c_th ];

    Rot = [ ...
        no_( 1 ) * no_( 1 ) * ( 1 - c_al ) + c_al, ...
        no_( 1 ) * no_( 2 ) * ( 1 - c_al ) - no_( 3 ) * s_al, ...
        no_( 1 ) * no_( 3 ) * ( 1 - c_al ) + no_( 2 ) * s_al; ...
        no_( 2 ) * no_( 1 ) * ( 1 - c_al ) + no_( 3 ) * s_al, ...
        no_( 2 ) * no_( 2 ) * ( 1 - c_al ) + c_al, ...
        no_( 2 ) * no_( 3 ) * ( 1 - c_al ) - no_( 1 ) * s_al; ...
        no_( 3 ) * no_( 1 ) * ( 1 - c_al ) - no_( 2 ) * s_al, ...
        no_( 3 ) * no_( 2 ) * ( 1 - c_al ) + no_( 1 ) * s_al, ...
        no_( 3 ) * no_( 3 ) * ( 1 - c_al ) + c_al ...
        ];    
    
    %% prepare time between two RF pulses
    
    % unique index
    
    lambda_time = 1;
    
    % constant gradient with 2 * pi dephasing per TR and dx_2pi
    
    p_scal = 1 ./ par.dx_2pi;
    
    % gradient moment in random direction
    
    p = Rot * [ p_scal; 0; 0 ];
    
    % gradient shape, as defined in the documentation (constant gradient)
    
    s = p;
    S = ( p * p' ) ./ 3;
    trace_S = sum( diag( S ) );
    
    %% initialize configuration model (idealized sequence)

    cm_0 = CoMoTk;      % without diffusion 
    cm = CoMoTk;        % with diffusion
    
    % mandatory tissue parameters
    
    cm_0.R1 = 1 / par.T1;
    cm_0.R2 = 1 / par.T2;
    cm_0.D = 0;
    
    cm.R1 = 1 / par.T1;
    cm.R2 = 1 / par.T2;
    
    if ( b_tensor )

        % calculate three orthonormal vectors e_j (j = 1,2,3)
        % with e_1 parallel to p
        
        e_1 = p ./ norm( p );
        
        tmp = null( e_1' );
        
        e_2 = tmp( :, 1 );
        e_3 = tmp( :, 2 );
        
        % generate diffusion tensor with one of the principal components
        % parallel to p (associated with par.D)
        % the diffusion coefficients along the other two principal
        % components (which should not matter) are chosen randomly
        
        D_orth = ( 0.5 + rand( 2, 1 ) ) .* par.D;
        
        cm.D = par.D .* ( e_1 * e_1' ) + D_orth( 1 ) .* ( e_2 * e_2' ) + D_orth( 2 ) .* ( e_3 * e_3' );
        
        % cm.D = par.D * eye( 3 );
    
    else
        
        cm.D = par.D;

    end
    
    % get default options
    
    options = cm.options;
    
    options.alloc_n = 1000;     % CoMoTk will allocate more, if needed
    options.alloc_d = 1;        % only one time interval
    options.epsilon = 0;
    
    % set new options
    
    cm_0.options = options;
    cm.options = options;

    % start with longitudinal magnetization
    
    cm_0.init_configuration ( [ 0; 0; 1 ] );
    cm.init_configuration ( [ 0; 0; 1 ] );
    
    %% approach steady state
    
    for i = 1 : num_TR
        
        if ( i > 1 )
        
            % time interval
            
            param = [];
            param.lambda = lambda_time;
            param.tau = par.TR;
            
            cm_0.time( param );
            
            param.p = p;
            
            if ( b_implicit )
                
                cm.time( param );
                
            else
                
                param.s = s;
                
                if ( b_tensor )
                    
                    param.S = S;
                    
                else
                    
                    param.S = trace_S;
                    
                end
                
                cm.time( param );
                
            end
            
        end
        
        % excitation pulse
        
        param = [];
        param.FlipAngle = fa_rad;
        param.Phase = 0;
                
        cm_0.RF( param );
        cm.RF( param );

    end
    
    %% allocate space for results
    
    m_xy_0 = zeros( 1, 2 * n_max + 1 );
    m_z_0 = zeros( 1, 2 * n_max + 1 );
    m_xy = zeros( 1, 2 * n_max + 1 );
    m_z = zeros( 1, 2 * n_max + 1 );

    % extract configurations for n = - n_max : n_max - 1
    
    for i = 1 : 2 * n_max + 1
        
        n = i - n_max - 1;
        
        b_n = cm_0.find( lambda_time, n );   % find location of configuration order n 
        
        m_xy_0( i ) = sqrt( 2 ) * cm_0.m( 1, b_n );
        m_z_0( i ) = cm_0.m( 2, b_n );
        
        m_xy( i ) = sqrt( 2 ) * cm.m( 1, b_n );
        m_z( i ) = cm.m( 2, b_n );
        
    end
    
    %% calculate analytical result
    
    freed_0 = diff_ssfp_freed( n_max, 3 * n_max, par.fa, par.TR, p_scal, par.T1, par.T2, 0 );
    freed = diff_ssfp_freed( n_max, 3 * n_max, par.fa, par.TR, p_scal, par.T1, par.T2, par.D );
    
    %% Show results
    
    n = - n_max : n_max;
    
    ax = subplot( 1, 1, 1 );
        
    semilogy( n, abs( m_xy ), '-b', n, abs( freed.xy ), '+b', n, abs( m_xy_0 ), '-r', n, abs( freed_0.xy ), '+r' );
    legend( [ 'CM  D = ', num2str( par.D ) ], [ 'Th.  D = ', num2str( par.D ) ], 'CM  D = 0', 'Th.  D = 0', 'Location', 'south' );
    title( 'transverse magnitude' );
    xlabel( 'configuration order n' );    
    
    xlim( [ (-par.n_max - 1) (par.n_max + 1) ] );
    ylim( [ 0.5 * min( abs( m_xy( : ) ) ) 2 * max( abs( freed_0.xy( : ) ) ) ] );
    
    width = 9;
    height = 9;
    
    set( gcf, 'Units', 'centimeters' );
    set( gcf, 'Position', [ 0, 0, width, height ] );
    set( gcf, 'Color', 'w' );
    
    delta = 0.01;
    
    ti = ax.TightInset;
    left = ti(1) + delta;
    bottom = ti(2) + delta;
    ax_width = 1 - ti(1) - ti(3) - 2 * delta;
    ax_height = 1 - ti(2) - ti(4) - 2 * delta;
    ax.Position = [left bottom ax_width ax_height ];
    
    %% show accuracy
    
    fprintf( 1, 'median( abs( CM / Freed - 1 ) )\n' );
    fprintf( 1, 'xy: %e\n', median( abs( m_xy / freed.xy - 1 ) ) );
    fprintf( 1, ' z: %e\n', median( abs( m_z / freed.z - 1 ) ) );
    fprintf( 1, 'max( abs( CM - Freed ) )\n' );
    fprintf( 1, 'xy: %e\n', max( abs( m_xy - freed.xy ) ) );
    fprintf( 1, ' z: %e\n', max( abs( m_z - freed.z ) ) );
    
end

function [ m_n,  n ] = diff_ssfp_freed( n_max, cf_max, alpha, TR, p, T1, T2, D )
% diff_ssfp_freed
%  
% Calculates transverse component of diffusion damped SSFP configuration vectors
% (immediately after RF pulse) according to Freed et al.
% (assumes constant gradient in each TR interval and RF phase = 0)
% calculated configuration orders = - n_max : n_max
%
% In:
% 
% n_max = maximal configuration order
% cf_max = depth of continued fraction
% alpha = flip angle [ deg ]
% TR = repetition time [ ms ]
% p = gradient moment [ 1 / um ]
% T1 = longitudinal relaxation time [ ms ]
% T2 = transverse relaxation time [ ms ]
% D = apparent diffusion coefficient [ um^2 / ms ]
%
% Out:
%
% m_n.xy = transverse component of configuration vector (normal convention)
% m_n.z = longitudinal component of configuration vector
% n = - n_max : n_max - 1

n = - n_max : n_max;
n_0 = n_max + 1;

cf = - cf_max : cf_max;
cf_0 = cf_max + 1;

E1 = exp( - TR / T1 - D .* p^2 .* TR .* cf.^2 );
E2 = exp( - TR / T2 - D .* p^2 .* TR .* ( cf.^2 + cf + 1 / 3 ) );

c_al = cosd( alpha );
s_al = sind( alpha );

A = 0.5 .* ( E1 - 1 ) .* ( 1 + c_al );
B = 0.5 .* ( E1 + 1 ) .* ( 1 - c_al );
C = E1 - c_al;

% the following quantities correspond to indices 1 : n_max - 1

n_n = - ...
      E2( cf_0 - 1 : -1 : 2 ) .* ...
      E2( cf_0 : end - 2 ) .* ...
      A( cf_0 + 1 : end - 1 ).^2 .* ...
      B( cf_0 : end - 2 ) ./ ...
      B( cf_0 + 1 : end - 1 );

e_n = - ...
      E2( cf_0 - 2 : -1 : 1 ) .* ...
      E2( cf_0 + 1 : end - 1 ) .* ...
      B( cf_0 + 1 : end - 1 ) .* ...
      C( cf_0 + 2 : end ) ./ ...
      B( cf_0 + 2 : end );

d_n = A( cf_0 + 1 : end - 1 ) - B( cf_0 + 1 : end - 1 ) - e_n;

% calculate the continued fraction

a_1 = 0;
a_2 = 1;
b_1 = 1;
b_2 = 0;

n_cut = length( d_n );

for i = 1 : n_cut  % n_max

    a = d_n( i ) * a_1 + n_n( i ) * a_2;
    b = d_n( i ) * b_1 + n_n( i ) * b_2;

    a_2 = a_1;
    a_1 = a;

    b_2 = b_1;
    b_1 = b;

end

x_1 = a / b;

% calculate r_1

r_1 = x_1 / ( E2( cf_0 - 1 ) *  B( cf_0 ) ) + ...
      E2( cf_0 ) *  C( cf_0 + 1 ) / B( cf_0 + 1 );

% allocate space for result

m_n.xy = zeros( 1, 2 * n_max );

% calculate m_xy_n for n == 0 (FID) and n == -1 (ECHO)

m_n.xy( n_0 ) = - s_al * ( 1 - E1( cf_0 ) ) ./ ...
    ( A( cf_0 ) - B( cf_0 ) + E2( cf_0 - 1 ) * C( cf_0 ) * r_1 );

m_n.xy( n_0 - 1 ) = - r_1 * m_n.xy( n_0 );

% calculate the other transverse configurations

for i = 1 : n_max
   
    m_n.xy( n_0 + i ) = ...
        ( E2( cf_0 + i - 1 ) * C( cf_0 + i ) * m_n.xy( n_0 + i - 1 ) + ...
          B( cf_0 + i ) * m_n.xy( n_0 - i ) ) ./ A( cf_0 + i );
    
      if ( i < n_max )
          
          m_n.xy( n_0 - i - 1 ) = ...
              ( - E2( cf_0 + i - 1 ) * B( cf_0 + i ) * m_n.xy( n_0 + i - 1 ) + ...
              ( A( cf_0 + i ) - B( cf_0 + i ) ) * m_n.xy( n_0 - i ) ) ./ ...
              ( E2( cf_0 - i - 1 ) * A( cf_0 + i ) );
          
      end
      
end

% RF phase == 0

m_n.z = zeros( 1, 2 * n_max + 1 );

m_n.z( n_0 ) = ( s_al * m_n.xy( n_0 ) - 1 + E1( cf_0 ) ) / C( cf_0 ); 

r_p = n_0 + 1 : n_0 + n_max;
r_n = n_0 - 1 : -1 : n_0 - n_max;
re_p = cf_0 + 1 : cf_0 + n_max;
re_n = cf_0 - 1 : -1 : cf_0 - n_max;

m_n.z( r_p ) = 0.5 .* s_al .* ( m_n.xy( r_p ) + m_n.xy( r_n ) ) ./ C( re_p );
m_n.z( r_n ) = 0.5 .* s_al .* ( m_n.xy( r_p ) + m_n.xy( r_n ) ) ./ C( re_n );

m_n.xy = - 1i .* m_n.xy;
m_n.z = real( m_n.z );

end
