%% GRE sequence
% This script compares an idealized versus realistic implementation of a GRE sequence

%% Simulation parameters
% (by default, all set by the user)

% Tissue properties

T1 = 1000; %input( 'T1 [ms] = \n' );
T2 = 100; %input( 'T2 [ms] = \n' );
D = 0; %input( 'D [um^2/ms] = \n' );         XXX diffusion NOT tested yet! XXX

% Sequence and measurement settings

num_TR = 25; %input( 'Number of TR cycles = \n' );

% TR 

TR = 5; %input( 'TR [ms] = \n' );

% B1

B1 = 1; %input( 'B1 = \n' );

%% RF pulse parameters
% to see the slice profile, check out "selective_excitation.m"
% (it is the same pulse)

% nominal flip angle (for B1 == 1 at center of the slice profile)

fa_deg = 30; %input( 'Excitation flip angle [deg] = \n' );

% phase increment (essentially off-resonance, used to shift bSSFP profile)

ph_inc_deg = 180; %input( 'Phase increment [deg] = \n' );

% phase difference increment (RF spoiling)

ph_diff_inc_deg = 0; %input( 'Phase difference increment [deg] = \n' );

% calculate phase cycling

ph_deg = zeros( num_TR, 1 );

for i = 2 : num_TR

  ph_deg( i ) = ph_deg( i - 1 ) + ph_inc_deg + ( i - 1 ) * ph_diff_inc_deg;
  
end

% pulse duration (has to be shorter than TR)

t_rf = TR;

while( t_rf >= TR )

    t_rf = 1; %input( 'RF pulse duration [ms] = \n' );
    
end

% slice thickness

sl_th_mm= 3; %input( 'slice thickness [mm] = \n' );

% number of RF support points 
%
% as in SLR optimization we approximate the RF pulse as an alternating train of
%
% *   n_tau + 1 instantaneous RF pulses, separated by
% *   n_tau     intervals of constant gradient amplitude

n_tau = 100; %input( 'number of RF pulse support intervals = \n' );

% the number of zero crossings (SINC pulse)
% n_zc == 1 means inclusion of the central peak only

n_zc = 5; %input( 'number of zero crossings per side = \n' );

% filtering on or off?

i = 1; %0;

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

cm_ideal = CoMoTk;

% mandatory tissue parameters

cm_ideal.R1 = 1 / T1;
cm_ideal.R2 = 1 / T2;
cm_ideal.D = D;

% further parameters

cm_ideal.B1 = B1;

% get default options

options = cm_ideal.options;

options.alloc_n = 1000;  % CoMoTk will allocate more, if needed

options.alloc_d = 1;     % 1 == half echo spacing

options.epsilon = 0;

% set new options

cm_ideal.options = options;

% start with longitudinal magnetization

cm_ideal.init_configuration ( [ 0; 0; 1 ] );

%% prepare time between end of RF pulse and echo

% unique index

mu_ideal_te = 1;

% duration

tau_ideal_te = 0.5 * TR;

%% initialize configuration model (real sequence)

cm_real = CoMoTk;

% mandatory tissue parameters

cm_real.R1 = 1 / T1;
cm_real.R2 = 1 / T2;
cm_real.D = D;

% further parameters

cm_real.B1 = B1;

% get default options

options = cm_real.options;

options.alloc_n = 10000;  % CoMoTk will allocate more, if needed

options.alloc_d = 2;     % 1 == RF pulse plateau
                                     % 2 == half echo spacing

options.epsilon = 1e-4;     % best accuracy

% set new options

cm_real.options = options;

% start with longitudinal magnetization

cm_real.init_configuration ( [ 0; 0; 1 ] );

%% prepare RF pulse plateau

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

%% prepare time between end of RF pulse and echo

% unique index

mu_real_te = 2;

% duration

tau_real_te = 0.5 * ( TR - t_rf );

% gradient moment == crusher and rephasing gradient

p_real_te = [ 0; 0; - p_sl / 2 ];

%% allocate space for results

m_ideal_SSFP = zeros( num_TR, 1 );
m_ideal_bSSFP = zeros( num_TR, 1 );
m_real_SSFP = zeros( num_TR, 1 );
m_real_bSSFP = zeros( num_TR, 1 );

%% (b)SSFP with instantaneous RF pulse

for i = 1 : num_TR

    % time from echo to next RF pulse
    % (reordering needed for the first readout)
    
    cm_ideal.time( mu_ideal_te, 'tau', tau_ideal_te );
    
    % excitation pulse
    
    cm_ideal.RF( fa_rad, ph_rad( i ) );
    
    % SSFP  : only the zero order configuration (= FID) contributes to the voxel signal,
    %         (assuming a crusher was present after the echo - even, if not simulated)
        
    b_n = cm_ideal.find( mu_ideal_te, 0 );
    
    % calculate the partial sum
    
    iso = cm_ideal.isochromat( 0, [], b_n );
        
    % save the echo
        
    m_ideal_SSFP( i ) = iso.xy;
        
    % time to TR / 2
        
    cm_ideal.time( mu_ideal_te, 'tau', tau_ideal_te );

    % bSSFP : all configurations contribute to the voxel signal
    
    iso = cm_ideal.isochromat( 0, [], [] );
        
    % save the echo
        
    m_ideal_bSSFP( i ) = iso.xy;
                
end

%% (b)SSFP with real RF pulse

for i = 1 : num_TR

    fprintf( 1, 'i = %d / %d\n', i, num_TR );
    
    % excitation pulse

    % first small pulse
    
    cm_real.RF( al_rad( 1 ), ph_rad( i ) );     
    
    for j = 1 : n_tau

        cm_real.time( mu_real_rf, 'tau', tau_real_rf, 'p', p_real_rf );
            
        % and the rest of the small pulses
        
        cm_real.RF( al_rad( j + 1 ), ph_rad( i ) );     
        
    end

    % SSFP  : only the zero order configuration (= FID) contributes to the voxel signal,
    %         (assuming a crusher was present after the echo - even, if not simulated)

    b_n = cm_real.find( mu_real_te, 0 );

    %         slice encoding direction requires a zero gradient moment
        
    b_n = b_n & reshape( cm_real.p_n( 3, : ) == 0, size( b_n ) );
    
    % calculate the partial sum
    
    iso = cm_real.isochromat( 0, [], b_n );
        
    % save the echo
        
    m_real_SSFP( i ) = iso.xy;
        
    % time to TR / 2
        
    cm_real.time( mu_real_te, 'tau', tau_real_te, 'p', p_real_te );

    % bSSFP : only the slice encoding direction requires a zero gradient moment
    
    b_n = cm_real.b_n & reshape( cm_real.p_n( 3, : ) == 0, size( cm_real.b_n ) );

    % calculate the partial sum
    
    iso = cm_real.isochromat( 0, [], b_n );
        
    % save the echo
        
    m_real_bSSFP( i ) = iso.xy;
        
    % time from echo to next RF pulse
        
    cm_real.time( mu_real_te, 'tau', tau_real_te, 'p', p_real_te );
        
end

%% Show results

te_SSFP = ( 0 : ( num_TR - 1 ) )' .* TR;
te_bSSFP = ( ( 1 : num_TR )' - 0.5 ) .* TR;

sc_SSFP = abs( m_ideal_SSFP( 1 ) / m_real_SSFP( 1 ) );
sc_bSSFP = abs( m_ideal_bSSFP( 1 ) / m_real_bSSFP( 1 ) );

subplot( 1, 2, 1 );
plot( te_SSFP, abs( m_ideal_SSFP ), te_SSFP, sc_SSFP .* abs( m_real_SSFP ) );
legend( 'ideal', 'real' );
title( 'SSFP' );

subplot( 1, 2, 2 );
plot( te_bSSFP, abs( m_ideal_bSSFP ), te_bSSFP, sc_bSSFP .* abs( m_real_bSSFP ) );
legend( 'ideal', 'real' );
title( 'bSSFP' );
