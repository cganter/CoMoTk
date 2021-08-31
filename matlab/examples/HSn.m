function out = HSn( in )
% HSn
%
% Returns hyperbolic secant (HSn) pulse in the linear rapid passage region
% based upon D. Idiyatullin et al., JMR 193 (2008) 267–273
% Supports gapped pulses with reduced duty cycle
%
% Input:
% 
% in.n = order of HSn pulse [no unit] (optional, default value == 1):
%
% At least *two* of the following three parameters (specifying duration and bandwidth)
% need to be specified.
% Their presence is checked in the following order
% (no check for inconsistent settings, if all three parameters are
% supplied):
%
% in.Tp = pulse duration [ms]
% in.bw = pulse bandwidth [kHz]
% in.R = time bandwidth product [no unit] (R = Tp * bw)
%
% flip angle for rapid passage HSn pulses:
%
% in.fa = flip angle [deg]
%
% the RF pulse is specified for rapid passage.
% Adiabatic pulses can be simulated, if the following optional value
% (default == 1)
%
% in.rapid 
%
% is set to a value < 1. (omega_1_max is divided by this value)
%
% Discretization of RF pulse:
%
% The RF pulse is assumed to be split into N_tot = N_seg * L_over parts
% each with constant B1+. 
% To satisfy the Nyquist condition, N_tot must be larger than R.
% For gapped pulses, the duty cycle d_c defines the fraction of B1+
% transmission and N_seg the number of gaps. 
% Consequently, HSn pulses without gap correspond to d_c == 1
% To avoid rounding errors, L_over * d_c should be an integer.
% We let gapped pulses start with a transmit part and end with a gap.
%
% in.L_over = oversampling factor [integer, no unit]
% in.N_seg = number of segments [no unit]
% in.d_c = duty cycle [no unit]
%
% Output:
%
% out.t = time points, related flip angles and phases [ms]
% out.alpha = flip angle train (length( alpha ) == L_over * N_seg + 1)
% out.phi = associated phases
% out.t_int = center times of small intervals
% out.omega_1_max = max. amplitude of omega_1
% out.omega_1 = time-dependent RF amplitude [rad * kHz]
% out.omega_RF = time-dependent RF frequency [rad * kHz]
% out.dw = dwell time (== Tp / N_seg)
% out.dt = oversampled time interval (== dw / L_over)
% out.n_B1 = number of transmit periods per dwell time (== L_over * d_c)
% out.n_gap = receive intervals per dwell time (== L_over - n_B1)
% out.a = frequency acceleration [kHz/ms] (== (d/dt omega_RF)/(2*pi))

% pulse order

if ( isfield( in, 'n' ) )

    n = in.n;
    
else
    
    n = 1;
    
end

% duration and bandwidth settings

if ( isfield( in, 'Tp' ) && isfield( in, 'bw' ) )
    
    Tp = in.Tp;
    bw = in.bw;
    % R = Tp * bw;
    
elseif ( isfield( in, 'Tp' ) && isfield( in, 'R' ) )
    
    Tp = in.Tp;
    R = in.R;
    bw = R / Tp;
    
elseif ( isfield( in, 'bw' ) && isfield( in, 'R' ) )
    
    bw = in.bw;
    R = in.R;
    Tp = R / bw;
    
else
    
    error( 'insufficient specification of duration and bandwidth.' );
    
end

% number of intervals into which the RF pulse is split

N_tot = in.N_seg * in.L_over;

% approximate flip angle in the rapid passage regime

fa_rad = in.fa * pi / 180;

% set dimensionless truncation factor

if ( isfield( in, 'beta' ) )

    beta = in.beta;
    
else
    
    beta = asech( 0.01 );
    
end

% location of gaps

out.n_B1 = round( in.d_c * in.L_over );     % number of B1+ transmit periods per segment
out.n_gap = in.L_over - out.n_B1;           % discretization of the gaps

b_gap = repmat( [ false( 1, out.n_B1 ), true( 1, out.n_gap ) ], [ 1, in.N_seg ] );

% time steps related to discretization

out.dw = Tp / in.N_seg;    % dwell time
out.dt = out.dw / in.L_over;   % discretization of the transmission periods

% Total time period of the pulse for the final hard pulse approximation
% At the t( j ), we place instantaneous rotations 
% with rotation angle out.alpha( j ) and phase out.phi( j ).

out.t = linspace( 0, Tp - out.dt, N_tot );

% center points of small intervals

out.t_int = out.t + 0.5 * out.dt;

% max. RF amplitude and (half) bandwidth [rad * kHz]

out.omega_1_max = fa_rad * sqrt( bw / Tp ) * beta .^ ( 0.5 / n );

if ( isfield( in, 'rapid' ) )
    
    out.omega_1_max = out.omega_1_max / in.rapid;
    
end

A = bw * pi;

% HSn driving function at center points

fn = sech( beta .* ( ( 2 / Tp ) .* out.t_int - 1 ) .^ n );

% time dependent RF amplitude

out.omega_1 = out.omega_1_max .* fn;

% which is set to zero in the gaps

out.omega_1( b_gap ) = 0;

% RF frequency

out.omega_RF = ( 2 * A ) .* ( cumsum( fn .^ 2 ) ./ sum( fn .^ 2 ) - 0.5 );

% frequeny acceleration

out.a = max( abs( out.omega_RF( 2 : end ) - out.omega_RF( 1 : end - 1 ) ) ) / ( 2 * pi * out.dt);

% flip angles in hard pulse approximation

out.alpha = out.omega_1 .* out.dt;

% compensate for reduced duty cycle

out.alpha = out.alpha ./ in.d_c;

% and the phases

out.phi = - out.dt .* cumsum( out.omega_RF );
        
end