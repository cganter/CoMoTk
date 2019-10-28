%% Compares diffusion of unbalanced SSFP with available analytical solution 
% derived by D. Freed et al. in: J Chem Phys (2001) 115, 4249-4258

%% Simulation parameters
% (by default, all set by the user)

par = [];
opt = [];
str = [];

par.T1 = 100;
par.T2 = 10;
par.D = 3;
par.TR = 5;
par.fa = 50;

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.TR = [];
opt.fa = [];

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.TR = '[ms] repetition time';
str.fa = '[deg] flip angle';

while ( true )
    
    [ par, sel ] = sfv( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    % number of TR cycles to approach steady state
    
    num_TR = 100 * ceil( par.T1 / par.TR );
    
    % convert to units, as expected by CoMoTk
    
    fa_rad = par.fa * pi / 180;
    
    % routine returns configuration orders n = - n_max : n_max - 1
    
    n_max = 10;
    
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
    
    %% initialize configuration model (idealized sequence)

    % isotropic version
    
    cm_scal_impl = CoMoTk;   % constant gradient assumed implicitly
    cm_scal_expl = CoMoTk;   % gradient shape provided explicitly
    
    % tensor version (setting tensor = D * eye( 3 ))
    
    cm_tens_impl = CoMoTk;   % constant gradient assumed implicitly
    cm_tens_expl = CoMoTk;   % gradient shape provided explicitly
        
    % mandatory tissue parameters
    
    cm_scal_impl.R1 = 1 / par.T1;
    cm_scal_impl.R2 = 1 / par.T2;
    cm_scal_impl.D = par.D;
    
    cm_scal_expl.R1 = 1 / par.T1;
    cm_scal_expl.R2 = 1 / par.T2;
    cm_scal_expl.D = par.D;
    
    cm_tens_impl.R1 = 1 / par.T1;
    cm_tens_impl.R2 = 1 / par.T2;
    cm_tens_impl.D = par.D * eye( 3 );
    
    cm_tens_expl.R1 = 1 / par.T1;
    cm_tens_expl.R2 = 1 / par.T2;
    cm_tens_expl.D = par.D * eye( 3 );
    
    % get default options
    
    options = cm_scal_impl.options;
    
    options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
    options.alloc_d = 1;        % only one time interval
    options.epsilon = 0;
    
    % set new options
    
    cm_scal_impl.options = options;
    cm_scal_expl.options = options;
    
    cm_tens_impl.options = options;
    cm_tens_expl.options = options;
    
    % start with longitudinal magnetization
    
    cm_scal_impl.init_configuration ( [ 0; 0; 1 ] );
    cm_scal_expl.init_configuration ( [ 0; 0; 1 ] );
    
    cm_tens_impl.init_configuration ( [ 0; 0; 1 ] );
    cm_tens_expl.init_configuration ( [ 0; 0; 1 ] );
    
    %% prepare time between two RF pulses
    
    % unique index
    
    mu_time = 1;
    
    % constant gradient with 2 * pi dephasing per TR and 100 um
    
    p_scal = 0.01;
    
    % gradient moment in random direction
    
    p = Rot * [ p_scal; 0; 0 ];
    
    % gradient shape, as defined in the documentation (constant gradient)
    
    s = p;
    S = ( p * p' ) ./ 3;
    trace_S = sum( diag( S ) );
    
    %% approach steady state
    
    for i = 1 : num_TR
        
        if ( i > 1 )
        
            % time interval -/+ providing the gradient shape explicitly
            % (implicitly, a constant gradient is assumed - so the result should be the same)
            
            param = [];
            param.mu = mu_time;
            param.tau = par.TR;
            param.p = p;
            
            cm_scal_impl.time( param );
            cm_tens_impl.time( param );

            param.s = s;
            param.S = trace_S;
            
            cm_scal_expl.time( param );

            param.S = S;
            
            cm_tens_expl.time( param );
        
        end
        
        % excitation pulse
               
        param = [];
        param.FlipAngle = fa_rad;
        param.Phase = 0;
                
        cm_scal_impl.RF( param );
        cm_scal_expl.RF( param );

        cm_tens_impl.RF( param );
        cm_tens_expl.RF( param );

    end
    
    %% allocate space for results
    
    m_xy_n_scal_impl = zeros( 1, 2 * n_max );
    m_xy_n_scal_expl = zeros( 1, 2 * n_max );

    m_xy_n_tens_impl = zeros( 1, 2 * n_max );
    m_xy_n_tens_expl = zeros( 1, 2 * n_max );
    
    % extract configurations for n = - n_max : n_max - 1
    
    for i = 1 : 2 * n_max
        
        n = i - n_max - 1;
        
        m_xy_n_scal_impl( i ) = sqrt( 2 ) * cm_scal_impl.m( 1, cm_scal_impl.find( mu_time, n ) );
        m_xy_n_scal_expl( i ) = sqrt( 2 ) * cm_scal_expl.m( 1, cm_scal_expl.find( mu_time, n ) );
        
        m_xy_n_tens_impl( i ) = sqrt( 2 ) * cm_tens_impl.m( 1, cm_tens_impl.find( mu_time, n ) );
        m_xy_n_tens_expl( i ) = sqrt( 2 ) * cm_tens_expl.m( 1, cm_tens_expl.find( mu_time, n ) );
        
    end
    
    %% calculate analytical result
    
    m_xy_n_freed = diff_ssfp_freed( n_max, 2 * n_max, par.fa, par.TR, p_scal, par.T1, par.T2, par.D );
    
    %% Show results
    
    n = - n_max : n_max - 1;
    
    subplot( 2, 2, 1 );
    
    plot( n, real( m_xy_n_scal_impl ), n, real( m_xy_n_scal_expl ), n, real( m_xy_n_freed ) );
    legend( 'implicit shape', 'explicit shape', 'analytical' );
    title( 'scalar, real part' );
    
    subplot( 2, 2, 3 );
    
    plot( n, imag( m_xy_n_scal_impl ), n, imag( m_xy_n_scal_expl ), n, imag( m_xy_n_freed ) );
    legend( 'implicit shape', 'explicit shape', 'analytical' );
    title( 'scalar, imaginary part' );
    
    subplot( 2, 2, 2 );
    
    plot( n, real( m_xy_n_tens_impl ), n, real( m_xy_n_tens_expl ), n, real( m_xy_n_freed ) );
    legend( 'implicit shape', 'explicit shape', 'analytical' );
    title( 'tensor, real part' );
    
    subplot( 2, 2, 4 );
    
    plot( n, imag( m_xy_n_tens_impl ), n, imag( m_xy_n_tens_expl ), n, imag( m_xy_n_freed ) );
    legend( 'implicit shape', 'explicit shape', 'analytical' );
    title( 'tensor, imaginary part' );
    
end

function [ m_xy_n, n ] = diff_ssfp_freed( n_max, cf_max, alpha, TR, p, T1, T2, D )
% diff_ssfp_freed
%  
% Calculates transverse component of diffusion damped SSFP configuration vectors
% (immediately after RF pulse) according to Freed et al.
% (assumes constant gradient in each TR interval and RF phase = 0)
% calculated configuration orders = - n_max : n_max - 1
%
% In:
% 
% n_max = maximal configuartion order
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
% m_xy_n = transverse component of configuration vector (normal convention)
% n = - n_max : n_max - 1

n = - n_max : n_max - 1;
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

for i = 1 : n_max - 1

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

m_xy_n = zeros( 1, 2 * n_max );

% calculate m_xy_n for n == 0 (FID) and n == -1 (ECHO)

m_xy_n( n_0 ) = - s_al * ( 1 - E1( cf_0 ) ) ./ ...
    ( A( cf_0 ) - B( cf_0 ) + E2( cf_0 - 1 ) * C( cf_0 ) * r_1 );

m_xy_n( n_0 - 1 ) = - r_1 * m_xy_n( n_0 );

% calculate the other transverse configurations

for i = 1 : n_max - 1
   
    m_xy_n( n_0 + i ) = ...
        ( E2( cf_0 + i - 1 ) * C( cf_0 + i ) * m_xy_n( n_0 + i - 1 ) + ...
          B( cf_0 + i ) * m_xy_n( n_0 - i ) ) ./ A( cf_0 + i );
    
    m_xy_n( n_0 - i - 1 ) = ...
        ( - E2( cf_0 + i - 1 ) * B( cf_0 + i ) * m_xy_n( n_0 + i - 1 ) + ...
          ( A( cf_0 + i ) - B( cf_0 + i ) ) * m_xy_n( n_0 - i ) ) ./ ...
        ( E2( cf_0 - i - 1 ) * A( cf_0 + i ) );
    
end

% RF phase == 0

m_xy_n = - 1i .* m_xy_n;

end
