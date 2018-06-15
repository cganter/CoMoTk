%% GRE sequence
% This script compares an idealized versus realistic implementation of a GRE sequence

%% Simulation parameters
% (by default, all set by the user)

% Tissue properties

T1 = 1000; %input( 'T1 [ms] = \n' );
T2 = 100; %input( 'T2 [ms] = \n' );
D = 0; %input( 'D [um^2/ms] = \n' );         XXX diffusion not tested yet! XXX

% Sequence and measurement settings

num_echos = 5; %input( 'Number of echos = \n' );
num_TR = 1; %input( 'Number of TR cycles = \n' );

echo_spacing = 10; %input( 'Echo spacing [ms] = \n' );

% TR (has to be larger than echo train)

TR = 0;

while( TR <= num_echos * echo_spacing )

    TR = 1000; %input( 'TR [ms] = \n' );
    
end

% to see diffusion effects due to 2\pi crushers

resolution_mm = 1; %input( 'In plane resolution [mm] = \n' );

B1 = 0.9; %input( 'B1 = \n' );

%% RF pulse parameters
% to see the slice profile, check out "selective_excitation.m"
% (it is the same pulse)

% nominal flip angle (for B1 == 1 at center of the slice profile)

fa_exc_deg = 90; %input( 'Excitation flip angle [deg] = \n' );
fa_ref_deg = 180; %input( 'Refocusing flip angle [deg] = \n' );

% corresponding phases
% CPMG condition == orthogonal rotation axes

ph_exc_deg = 0; %input( 'Excitation phase [deg] = \n' );
ph_ref_deg = 90; %input( 'Refocusing phase [deg] = \n' );

% pulse duration (has to be shorter than echo spacing)

t_rf = echo_spacing;

while( t_rf >= echo_spacing )

    t_rf = 2; %input( 'RF pulse duration [ms] = \n' );
    
end

% slice thickness

sl_th_mm= 3; %input( 'slice thickness [mm] = \n' );

% number of RF support points 
%
% as in SLR optimization we approximate the RF pulse as an alternating train of
%
% *   n_tau + 1 instantaneous RF pulses, separated by
% *   n_tau     intervals of constant gradient amplitude

n_tau = 300; %input( 'number of RF pulse support intervals = \n' );

% the number of zero crossings (SINC pulse)
% n_zc == 1 means inclusion of the central peak only

n_zc = 10; %input( 'number of zero crossings per side = \n' );

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

fa_exc_rad = fa_exc_deg * pi / 180;
fa_ref_rad = fa_ref_deg * pi / 180;
ph_exc_rad = ph_exc_deg * pi / 180;
ph_ref_rad = ph_ref_deg * pi / 180;
sl_th_um = 1000 * sl_th_mm;
resolution_um = 1000 * resolution_mm;

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

al_exc_rad = al .* ( fa_exc_rad / sum( al ) );
al_ref_rad = al .* ( fa_ref_rad / sum( al ) );

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

options.alloc_d = 2;     % 1 == half echo spacing with crusher
                         % 2 == time between last echo readout and next excitation pulse 
                         %      (includes ideal spoiler)

options.epsilon = 0;

% set new options

cm_ideal.options = options;

% start with longitudinal magnetization

cm_ideal.init_configuration ( [ 0; 0; 1 ] );

%% prepare time between end of RF pulse and echo

% unique index

mu_ideal_te = 1;

% duration

tau_ideal_te = 0.5 * echo_spacing;

% gradient moment == crusher

p_ideal_te = [ 2 * pi / resolution_um; 0; 0 ];

%% prepare time between last echo readout and next excitation pulse

% unique index

mu_ideal_tr = 2;

% duration

tau_ideal_tr = TR - num_echos * echo_spacing;

% gradient moment == Inf (ideal spoiler, even if D == 0)

p_ideal_tr = [ Inf; 0; 0 ];

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

options.alloc_d = 3;     % 1 == RF pulse plateau
                         % 2 == time between end of RF pulse and echo
                         % 3 == time between last echo readout and next excitation pulse 
                         %      (includes ideal spoiler)

% options.verbose = true;                         
                         
options.epsilon = 0;     % best accuracy

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

tau_real_te = 0.5 * ( echo_spacing - t_rf );

% gradient moment == crusher and rephasing gradient

p_real_te = [ 2 * pi / resolution_um; 0; - p_sl / 2 ];

%% in the first interval, we need two slice rephasing gradients

% unique index

mu_real_first_te = 4;

% duration

tau_real_first_te = 0.5 * ( echo_spacing - t_rf );

% gradient moment == crusher and 2 * rephasing gradient

p_real_first_te = [ 2 * pi / resolution_um; 0; - p_sl ];

%% prepare time between last echo readout and next excitation pulse

% unique index

mu_real_tr = 3;

% duration

tau_real_tr = TR - num_echos * echo_spacing;

% gradient moment == Inf (ideal spoiler, even if D == 0)

p_real_tr = [ Inf; 0; 0 ];

%% allocate space for results

m_ideal = zeros( num_echos, num_TR );
m_real = zeros( num_echos, num_TR );

%% TSE with instantaneous RF pulse

for i = 1 : num_TR

    % excitation pulse
    
    cm_ideal.RF( fa_exc_rad, ph_exc_rad );

    for j = 1 : num_echos       
        
        % time to refucusing pulse
        
        cm_ideal.time( mu_ideal_te, 'tau', tau_ideal_te, 'p', p_ideal_te );
        
        % refocusing pulse
        
        cm_ideal.RF( fa_ref_rad, ph_ref_rad );
        
        % time to echo
        
        cm_ideal.time( mu_ideal_te, 'tau', tau_ideal_te, 'p', p_ideal_te );
    
        % select the correct configuration:
        % only the zero order configuration along the crusher direction contributes to the voxel signal
        
        b_n = cm_ideal.find( mu_ideal_te, 0 );
        
        % calculate the partial sum
        
        iso = cm_ideal.isochromat( 0, [], b_n );
        
        % save the echo
        
        m_ideal( j, i ) = iso.xy;
        
    end
        
    if ( i < num_TR )

        % time to next excitation pulse
                
        cm_ideal.time( mu_ideal_tr, 'tau', tau_ideal_tr, 'p', p_ideal_tr );
    
    end
        
end

%% TSE with real RF pulse

for i = 1 : num_TR

    % excitation pulse

    % first small pulse
    
    cm_real.RF( al_exc_rad( 1 ), ph_exc_rad );     
    
    for k = 1 : n_tau

        cm_real.time( mu_real_rf, 'tau', tau_real_rf, 'p', p_real_rf );
            
        % and the rest of the small pulses
        
        cm_real.RF( al_exc_rad( k + 1 ), ph_exc_rad );     
        
    end

    for j = 1 : num_echos
        
        fprintf( 'i = %d, j = %d\n', i, j );
        
        % time to refucusing pulse
        
        if ( j == 1 )
        
            cm_real.time( mu_real_first_te, 'tau', tau_real_first_te, 'p', p_real_first_te );
    
        else
            
            cm_real.time( mu_real_te, 'tau', tau_real_te, 'p', p_real_te );

        end
        
        % refocusing pulse
        
        % first small pulse
        
        cm_real.RF( al_ref_rad( 1 ), ph_ref_rad );     
        
        for k = 1 : n_tau

            % small time interval
            
            cm_real.time( mu_real_rf, 'tau', tau_real_rf, 'p', p_real_rf );
                
            % small RF pulse
            
            cm_real.RF( al_ref_rad( k + 1 ), ph_ref_rad );     
            
        end

        % time to echo
        
        cm_real.time( mu_real_te, 'tau', tau_real_te, 'p', p_real_te );
    
        % select the correct configuration:
        % only the zero order configuration along the crusher direction contributes to the voxel signal
        
        % b_n = cm_real.find( mu_real_te, 0 );
        
        % additionally the gradient moments along the slice normal must vanish
        
        %         b_n = cm_real.b_n & reshape( cm_real.p_n( 1, : ) == 0, size( cm_real.b_n ) );
        %         b_n = b_n & reshape( cm_real.p_n( 3, : ) == 0, size( b_n ) );

        b_n = cm_real.b_n & reshape( ~any( cm_real.p_n ), size( cm_real.b_n ) );
        
        % calculate the partial sum
        
        iso = cm_real.isochromat( 0, [], b_n );
        
        % save the echo
        
        m_real( j, i ) = iso.xy;
                
    end

    if ( i < num_TR )
    
        % time to next excitation pulse
                
        cm_real.time( mu_real_tr, 'tau', tau_real_tr, 'p', p_real_tr );
        
    end
        
end

%% Compare the decays of the last cycle

te = ( 1 : num_echos )' .* echo_spacing;

sc = abs( m_ideal( 1, end ) / m_real( 1, end ) );

subplot( 1, 2, 1 );
plot( te, abs( m_ideal( :, end ) ), te, sc .* abs( m_real( :, end ) ) );
legend( 'ideal TSE', 'real TSE' );

subplot( 1, 2, 2 );
semilogy( te, abs( m_ideal( :, end ) ), te, sc .* abs( m_real( :, end ) ) );
legend( 'ideal TSE', 'real TSE' );
