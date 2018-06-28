%% Compares diffusion of unbalanced SSFP with available analytical solution 
% derived by D. Freed et al. in: J Chem Phys (2001) 115, 4249-4258

%% Simulation parameters
% (by default, all set by the user)

% Tissue properties

T1 = input( 'T1 [ms] = \n' );
T2 = input( 'T2 [ms] = \n' );
D = input( 'D [um^2/ms] = \n' );

% Sequence and measurement settings

% TR 

TR = input( 'TR [ms] = \n' );

% number of TR cycles to approach steady state

num_TR = 10 * ceil( T1 / TR );

% flip angle

fa_deg = input( 'Excitation flip angle [deg] = \n' );

% convert to units, as expected by CoMoTk

fa_rad = fa_deg * pi / 180;

% routine returns configuration orders n = - n_max : n_max - 1

n_max = 10;

%% initialize configuration model (idealized sequence)

cm_impl = CoMoTk;   % constant gradient assumed implicitly
cm_expl = CoMoTk;   % gradient shape provided explicitly

% mandatory tissue parameters

cm_impl.R1 = 1 / T1;
cm_impl.R2 = 1 / T2;
cm_impl.D = D;

cm_expl.R1 = 1 / T1;
cm_expl.R2 = 1 / T2;
cm_expl.D = D;

% get default options

options = cm_impl.options;

options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
options.alloc_d = 1;        % only one time interval
options.epsilon = 0;        % do not discard any modes

% set new options

cm_impl.options = options;
cm_expl.options = options;

% start with longitudinal magnetization

cm_impl.init_configuration ( [ 0; 0; 1 ] );
cm_expl.init_configuration ( [ 0; 0; 1 ] );

%% prepare time between two RF pulses

% unique index

mu_time = 1;

% constant gradient with 2 * pi dephasing per TR and 100 um

p_scalar = 0.01;

% gradient moment

p = [ p_scalar; 0; 0 ];

% gradient shape, as defined in the documentation (constant gradient)

s = [ p_scalar; 0; 0; p_scalar^2 / 3 ];

%% approach steady state

for i = 1 : num_TR

    % time interval -/+ providing the gradient shape explicitly
    % (implicitly, a constant gradient is assumed - so the result should be the same)
        
    cm_impl.time( mu_time, 'tau', TR, 'p', p );
    cm_expl.time( mu_time, 'tau', TR, 'p', p, 's', s );

    % excitation pulse
    
    cm_impl.RF( fa_rad, 0 );
    cm_expl.RF( fa_rad, 0 );
    
end

%% allocate space for results

m_xy_n_impl = zeros( 1, 2 * n_max );
m_xy_n_expl = zeros( 1, 2 * n_max );

% extract configurations for n = - n_max : n_max - 1

for i = 1 : 2 * n_max
   
    n = i - n_max - 1;

    m_xy_n_impl( i ) = sqrt( 2 ) * cm_impl.m( 1, cm_impl.find( mu_time, n ) );
    m_xy_n_expl( i ) = sqrt( 2 ) * cm_expl.m( 1, cm_expl.find( mu_time, n ) );
    
end

%% calculate analytical result

m_xy_n_freed = diff_ssfp_freed( n_max, 2 * n_max, fa_deg, TR, p_scalar, T1, T2, D );

%% Show results

n = - n_max : n_max - 1;

subplot( 1, 2, 1 );

plot( n, real( m_xy_n_impl ), n, real( m_xy_n_expl ), n, real( m_xy_n_freed ) );
legend( 'implicit shape', 'explicit shape', 'analytical' );
title( 'real part' );

subplot( 1, 2, 2 );

plot( n, imag( m_xy_n_impl ), n, imag( m_xy_n_expl ), n, imag( m_xy_n_freed ) );
legend( 'implicit shape', 'explicit shape', 'analytical' );
title( 'imaginary part' );

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
