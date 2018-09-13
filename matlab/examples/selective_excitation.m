%% Simple selective excitation
% This script shows how real RF pulses can be implemented.
% The example implements a simple amplitude modulated pulse without relaxation effects.
% The generalization to phase modulated and/or adiabatic pulses works along the same lines.

%% Tissue: we don't include relaxation or diffusion effects
% (if needed, just set the tissue parameters accordingly)

R1 = input( 'R1 [1/ms] = \n' );
R2 = input( 'R2 [1/ms] = \n' );
D = input( 'D [um^2/ms] = \n' );

%% RF pulse parameters
% flip angle (at center of the slice profile)

alpha_deg = input( 'flip angle [deg] = \n' );

% pulse duration (has to be shorter than echo spacing)

t_rf = input( 'RF pulse duration [ms] = \n' );
    
% slice thickness

sl_th_mm = input( 'slice thickness [mm] = \n' );

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

alpha_rad = alpha_deg * pi / 180;
sl_th_um = 1000 * sl_th_mm;

%% now we set the pulse profile 
% as described in the Handbook of MRI for a SINC pulse +/- Hamming filter

% duration of small intervals

tau_rf = t_rf / n_tau;

% support points and flip angles of instantaneous RF pulses (\propto B1 profile)

t = n_zc .* linspace( -1, 1, n_tau + 1 );

if ( hamming_filter )

    alpha = ( 0.54 + 0.46 .* cos( pi .* t ./ n_zc ) ) .* sinc( t );
    
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
x( 3, : ) = linspace( - sl_th_mm, sl_th_mm, n_sl );

%% initialize configuration model

cm = CoMoTk;

% mandatory tissue parameters

cm.R1 = R1;
cm.R2 = R2;
cm.D = D;

% get default options

options = cm.options;

options.alloc_n = 1000;
options.alloc_d = 2;     % == slice selection and rephasing gradient
options.epsilon = 0;     % best accuracy

% set new options

cm.options = options;

% start with longitudinal magnetization

cm.init_configuration ( [ 0; 0; 1 ] );

%% calculate slice selection and rephasing gradient moment

% RF pulse bandwidth

f = 2 * n_zc / t_rf;

% total moment of slice selection gradient

p_sl = 2 * pi * f * t_rf / sl_th_mm;

% split into n_tau intervals

p_rf = [ 0; 0; p_sl / n_tau ];

% rephasing gradient moment (for simplicity taken as - 0.5 * p_sl)
% as discussed in the Handbook of MRI (and shown by the results),
% this choice is reasonable but not optimal.

p_refoc = [ 0; 0; - p_sl / 2 ];

% assign unique handles for the two time intervals

mu_rf = 1;
mu_refoc = 2;

%% execute the RF pulse

for i = 1 : n_tau

    if ( i == 1 )

        % we start with the first small RF pulse
        
        cm.RF( alpha( i ), phase_rad ); 
        
        % in the first time interval, we specify the parameters
        
        cm.time( mu_rf, 'tau', tau_rf, 'p', p_rf );
    
    else
    
        % in subsequent calls, we only need the handle
        
        cm.time( mu_rf );
    
    end
    
    % and the rest of the small pulses
    
    cm.RF( alpha( i + 1 ), phase_rad );
    
end

%% collect the slice profile at locations x
% prior to the rephasing gradient, the phase will vary across the slice

m_iso_rf = zeros( 3, n_sl );

for i = 1 : n_sl
    
    % we sum over all configurations (third argument == [])

    iso = cm.isochromat( omega, x( :, i ), [] );
    
    m_iso_rf( 1, i ) = real( iso.xy );
    m_iso_rf( 2, i ) = imag( iso.xy );
    m_iso_rf( 3, i ) = real( iso.z );

end

%% apply the rephasing gradient
% in absense of relaxation, the duration does not matter
    
tau_refoc = 1; 

cm.time( mu_refoc, 'tau', tau_refoc, 'p', p_refoc );

%% collect the slice profile at locations x
% now the phase should be essentially constant across the slice

m_iso_refoc = zeros( 3, n_sl );

for i = 1 : n_sl
    
    % a restriction of configurations is still not neccessary

    iso = cm.isochromat( omega, x( :, i ), [] );
    
    m_iso_refoc( 1, i ) = real( iso.xy );
    m_iso_refoc( 2, i ) = imag( iso.xy );
    m_iso_refoc( 3, i ) = real( iso.z );

end

%% look at the results

subplot( 1, 3, 1);
plot( t_rf .* linspace( -0.5, 0.5, n_tau + 1 ), alpha ./ max( alpha ) );
xlim( [ - 0.5 * t_rf 0.5 * t_rf ] );
title( 'B^+_1(t)' );

subplot( 1, 3, 2 );
plot( x( 3, : ), m_iso_rf( 1, : ), x( 3, : ), m_iso_rf( 2, : ), x( 3, : ), m_iso_rf( 3, : ) );
legend( 'm_x', 'm_y', 'm_z' );
title( 'after RF pulse' );

subplot( 1, 3, 3 );
plot( x( 3, : ), m_iso_refoc( 1, : ), x( 3, : ), m_iso_refoc( 2, : ), x( 3, : ), m_iso_refoc( 3, : ) );
legend( 'm_x', 'm_y', 'm_z' );
title( 'after rephasing gradient' );
