%% GRE sequence
% This script compares an idealized versus realistic implementation of a GRE sequence

%% Simulation parameters
% (by default, all set by the user)

% Tissue properties

T1 = input( 'T1 [ms] = \n' );
T2 = input( 'T2 [ms] = \n' );
D = input( 'D [um^2/ms] = \n' );

% Sequence and measurement settings

num_TR = input( 'Number of TR cycles = \n' );

% TR 

TR = input( 'TR [ms] = \n' );

% pulse duration (has to be shorter than TR)

t_rf = TR;

while( t_rf <= 0 || t_rf >= TR )

    t_rf = input( 'RF pulse duration [ms] = \n' );
    
end

% TE relevant for unbalanced SSFP, for bSSFP TE = TR/2

TE = 0;

while( TE < t_rf / 2 || TE > TR - t_rf / 2 )

    TE = input( 'TE for unbalanced SSFP [ms] = \n' );
    
end

% B1

B1 = input( 'B1 = \n' );

%% RF pulse parameters
% to see the slice profile, check out "selective_excitation.m"
% (it is the same pulse)

% nominal flip angle (for B1 == 1 at center of the slice profile)

fa_deg = input( 'Excitation flip angle [deg] = \n' );

% phase increment (essentially off-resonance, used to shift bSSFP profile)

ph_inc_deg = input( 'Phase increment [deg] = \n' );

% phase difference increment (RF spoiling)

ph_diff_inc_deg = input( 'Phase difference increment [deg] = \n' );

% calculate phase cycling

ph_deg = zeros( num_TR, 1 );

for i = 2 : num_TR

  ph_deg( i ) = ph_deg( i - 1 ) + ph_inc_deg + ( i - 1 ) * ph_diff_inc_deg;
  
end

% slice thickness

sl_th_mm= input( 'slice thickness [mm] = \n' );

% number of RF support points 
%
% as in SLR optimization we approximate the RF pulse as an alternating train of
%
% *   n_tau + 1 instantaneous RF pulses, separated by
% *   n_tau     intervals of constant gradient amplitude

n_tau = input( 'number of RF pulse support intervals = \n' );

% the number of zero crossings (SINC pulse)
% n_zc == 1 means inclusion of the central peak only

n_zc = input( 'number of zero crossings per side = \n' );

% filtering on or off?

i = 0;

while ( i ~= 1 && i ~= 2 )

    fprintf( 1, 'Hamming filter?\n' );
    fprintf( 1, '1 == yes, 2 == no\n' );
    i = input( 'i = \n' );
    
end
 
if ( i == 1 )
    
    hamming_filter = true;
    
else
    
    hamming_filter = false;
    
end

% convert to units, as expected by CoMoTk

fa_rad = fa_deg * pi / 180;
ph_rad = ph_deg * pi / 180;
sl_th_um = 1000 * sl_th_mm;

%% set the RF pulse profile 
% as described in the Handbook of MRI for a SINC pulse +/- Hamming filter

% Number of RF support points
%
% as in SLR optimization we approximate the RF pulse as an alternating train of
%
% *   n_tau + 1 instantaneous RF pulses, separated by
% *   n_tau     intervals of constant gradient amplitude

% duration of small intervals

tau_rf = t_rf / n_tau;

% support points and flip angles of instantaneous RF pulses (\propto B1 profile)

t = n_zc .* linspace( -1, 1, n_tau + 1 );

if ( hamming_filter )

    al = ( 0.54 + 0.46 .* cos( pi .* t ./ n_zc ) ) .* sinc( t );
    
else
    
    al = sinc( t );
    
end

% normalization: at the center of the slice profile, we need the desired flip angle

al_rad = al .* ( fa_rad / sum( al ) );

%% initialize configuration model (idealized sequence)

% ======== unbalanced ========

cm_SSFP_ideal = CoMoTk;

% mandatory tissue parameters

cm_SSFP_ideal.R1 = 1 / T1;
cm_SSFP_ideal.R2 = 1 / T2;
cm_SSFP_ideal.D = D;

% further parameters

cm_SSFP_ideal.B1 = B1;

% set options

options = cm_SSFP_ideal.options;
options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
options.alloc_d = 2;        % before and after echo
options.epsilon = 0;
cm_SSFP_ideal.options = options;

% start with longitudinal magnetization

cm_SSFP_ideal.init_configuration ( [ 0; 0; 1 ] );

% prepare time between end of RF pulse and echo
% unique index

mu_SSFP_ideal_pre = 1;
mu_SSFP_ideal_post = 2;

% duration

tau_SSFP_ideal_pre = TE;
tau_SSFP_ideal_post = TR - TE;

% ======== balanced ========

cm_bSSFP_ideal = CoMoTk;

% mandatory tissue parameters

cm_bSSFP_ideal.R1 = 1 / T1;
cm_bSSFP_ideal.R2 = 1 / T2;
cm_bSSFP_ideal.D = D;

% further parameters

cm_bSSFP_ideal.B1 = B1;

% set options

options = cm_bSSFP_ideal.options;
options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
options.alloc_d = 1;        % TE = TR / 2
options.epsilon = 0;
cm_bSSFP_ideal.options = options;

% start with longitudinal magnetization

cm_bSSFP_ideal.init_configuration ( [ 0; 0; 1 ] );

% prepare time between end of RF pulse and echo
% unique index

mu_bSSFP_ideal = 1;

% duration

tau_bSSFP_ideal = TR / 2;

%% initialize configuration model (real sequence)

% ======== prepare RF pulse plateau ========

% unique index

mu_real_rf = 1;

% duration

tau_real_rf = t_rf / n_tau;

% RF pulse bandwidth

f = 2 * n_zc / t_rf;

% total moment of slice selection gradient

p_sl = 2 * pi * f * t_rf / sl_th_um;

% split into n_tau intervals

p_real_rf = [ 0; 0; p_sl / n_tau ];

% ======== unbalanced ========

cm_SSFP_real = CoMoTk;

% mandatory tissue parameters

cm_SSFP_real.R1 = 1 / T1;
cm_SSFP_real.R2 = 1 / T2;
cm_SSFP_real.D = D;

% further parameters

cm_SSFP_real.B1 = B1;

% set options

options = cm_SSFP_real.options;
options.alloc_n = 10000;  % CoMoTk will allocate more, if needed
options.alloc_d = 3;          % 1: RF pulse plateau
                                          % 2,3: before and after echo
options.epsilon = 0;
cm_SSFP_real.options = options;

% start with longitudinal magnetization

cm_SSFP_real.init_configuration ( [ 0; 0; 1 ] );

% prepare times between RF pulse and echo
% unique index

mu_SSFP_real_pre = 2;
mu_SSFP_real_post = 3;

% duration

tau_SSFP_real_pre = TE - t_rf / 2;
tau_SSFP_real_post = TR - TE - t_rf / 2;

% gradient moment == crusher and rephasing gradient

p_SSFP_real_pre = [ 0; 0; - p_sl / 2 ];
p_SSFP_real_post = [ 0; 0; - p_sl / 2 ];

% ======== balanced ========

cm_bSSFP_real = CoMoTk;

% mandatory tissue parameters

cm_bSSFP_real.R1 = 1 / T1;
cm_bSSFP_real.R2 = 1 / T2;
cm_bSSFP_real.D = D;

% further parameters

cm_bSSFP_real.B1 = B1;

% set options

options = cm_bSSFP_real.options;
options.alloc_n = 10000;  % CoMoTk will allocate more, if needed
options.alloc_d = 2;          % 1: RF pulse plateau
                                          % 2: TE = TR / 2
options.epsilon = 0;
cm_bSSFP_real.options = options;

% start with longitudinal magnetization

cm_bSSFP_real.init_configuration ( [ 0; 0; 1 ] );

% prepare times between RF pulse and echo
% unique index

mu_bSSFP_real = 2;

% duration

tau_bSSFP_real = ( TR - t_rf ) / 2;

% gradient moment == crusher and rephasing gradient

p_bSSFP_real = [ 0; 0; - p_sl / 2 ];

%% allocate space for results

m_SSFP_ideal = zeros( num_TR, 1 );
m_bSSFP_ideal = zeros( num_TR, 1 );
m_SSFP_real = zeros( num_TR, 1 );
m_bSSFP_real = zeros( num_TR, 1 );

for i = 1 : num_TR

    fprintf( 1, 'i = %d / %d\n', i, num_TR );
    
    %% (b)SSFP with instantaneous RF pulse

    % ======== unbalanced ========
    
    % time after echo
    % placed here, to initialize the crusher time interval
    % otherwise b_n cannot bet set for i == 1 below
    % since the initial state is in equilibrium for i == 1, the magnetization is not changed
        
    cm_SSFP_ideal.time( mu_SSFP_ideal_post, 'tau', tau_SSFP_ideal_post );
        
    % excitation pulse
    
    cm_SSFP_ideal.RF( fa_rad, ph_rad( i ) );
        
    % time to echo
        
    cm_SSFP_ideal.time( mu_SSFP_ideal_pre, 'tau', tau_SSFP_ideal_pre );
    
    % SSFP  : only the zero order configuration (= FID) contributes to the voxel signal,
    %         (assuming a crusher was present after the echo - even, if not simulated)
        
    b_n = cm_SSFP_ideal.find( mu_SSFP_ideal_post, 0 );
    
    % calculate the partial sum
    
    iso = cm_SSFP_ideal.isochromat( 0, [], b_n );
        
    % save the echo
        
    m_SSFP_ideal( i ) = iso.xy;
 
    % ======== balanced ========
    
    % time after echo
    
    cm_bSSFP_ideal.time( mu_bSSFP_ideal, 'tau', tau_bSSFP_ideal );
    
    % excitation pulse
    
    cm_bSSFP_ideal.RF( fa_rad, ph_rad( i ) );
    
    % time to echo
        
    cm_bSSFP_ideal.time( mu_bSSFP_ideal, 'tau', tau_bSSFP_ideal );

    % bSSFP : all configurations contribute to the voxel signal
    
    iso = cm_bSSFP_ideal.isochromat( 0, [], [] );
        
    % save the echo
        
    m_bSSFP_ideal( i ) = iso.xy;
                
    %% (b)SSFP with real RF pulse

    % ======== unbalanced ========
    
    % time after echo
    
    cm_SSFP_real.time( mu_SSFP_real_post, 'tau', tau_SSFP_real_post, 'p', p_SSFP_real_post );
    
    % excitation pulse

    % first small pulse
    
    cm_SSFP_real.RF( al_rad( 1 ), ph_rad( i ) );     
    
    for j = 1 : n_tau

        cm_SSFP_real.time( mu_real_rf, 'tau', tau_real_rf, 'p', p_real_rf );
            
        % and the rest of the small pulses
        
        cm_SSFP_real.RF( al_rad( j + 1 ), ph_rad( i ) );     
        
    end

    % time to echo
        
    cm_SSFP_real.time( mu_SSFP_real_pre, 'tau', tau_SSFP_real_pre, 'p', p_SSFP_real_pre );
    
    % SSFP  : only the zero order configuration (= FID) contributes to the voxel signal,
    %         (assuming a crusher was present after the echo - even, if not simulated)

    b_n = cm_SSFP_real.find( mu_SSFP_real_post, 0 );

    % slice encoding direction requires a zero gradient moment
        
    b_n = b_n & reshape( cm_SSFP_real.p_n( 3, : ) == 0, size( b_n ) );
    
    % calculate the partial sum
    
    iso = cm_SSFP_real.isochromat( 0, [], b_n );
        
    % save the echo
        
    m_SSFP_real( i ) = iso.xy;
        
    % ======== balanced ========
    
    % time after echo
    
    cm_bSSFP_real.time( mu_bSSFP_real, 'tau', tau_bSSFP_real, 'p', p_bSSFP_real );
    
    % excitation pulse

    % first small pulse
    
    cm_bSSFP_real.RF( al_rad( 1 ), ph_rad( i ) );     
    
    for j = 1 : n_tau

        cm_bSSFP_real.time( mu_real_rf, 'tau', tau_real_rf, 'p', p_real_rf );
            
        % and the rest of the small pulses
        
        cm_bSSFP_real.RF( al_rad( j + 1 ), ph_rad( i ) );     
        
    end

    % time to echo
        
    cm_bSSFP_real.time( mu_bSSFP_real, 'tau', tau_bSSFP_real, 'p', p_bSSFP_real );
    
    % bSSFP : only the slice encoding direction requires a zero gradient moment
    
    b_n = cm_bSSFP_real.b_n & reshape( cm_bSSFP_real.p_n( 3, : ) == 0, size( cm_bSSFP_real.b_n ) );

    % calculate the partial sum
    
    iso = cm_bSSFP_real.isochromat( 0, [], b_n );
        
    % save the echo
        
    m_bSSFP_real( i ) = iso.xy;
                
end

%% Show results

te_SSFP = TE + ( 0 : num_TR - 1 )' .* TR;
te_bSSFP = TR / 2 + ( 0 : num_TR - 1 )' .* TR;

sc_SSFP = abs( m_SSFP_ideal( 1 ) / m_SSFP_real( 1 ) );
sc_bSSFP = abs( m_bSSFP_ideal( 1 ) / m_bSSFP_real( 1 ) );

subplot( 1, 2, 1 );
plot( te_SSFP, abs( m_SSFP_ideal ), te_SSFP, sc_SSFP .* abs( m_SSFP_real ) );
legend( 'ideal', 'real' );
title( 'SSFP' );

subplot( 1, 2, 2 );
plot( te_bSSFP, abs( m_bSSFP_ideal ), te_bSSFP, sc_bSSFP .* abs( m_bSSFP_real ) );
legend( 'ideal', 'real' );
title( 'bSSFP' );
