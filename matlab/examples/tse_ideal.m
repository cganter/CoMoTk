%% TSE sequence
% This script compares idealized TSE with and without CPMG condition

%% Simulation parameters
% (by default, all set by the user)

% Tissue properties

T1 = 1000; %input( 'T1 [ms] = \n' );
T2 = 100; %input( 'T2 [ms] = \n' );
D = 0; %input( 'D [um^2/ms] = \n' );

% Sequence and measurement settings

num_echos = 25; %input( 'Number of echos = \n' );
num_cycles = 1; %input( 'Number of TR cycles = \n' );

echo_spacing = 10; %input( 'Echo spacing [ms] = \n' );

% TR (has to be larger than echo train)

TR = 0;

while( TR <= num_echos * echo_spacing )

    TR = 1000; %input( 'TR [ms] = \n' );
    
end

% to see diffusion effects due to 2 * pi crushers

resolution_mm = 1; %input( 'In plane resolution [mm] = \n' );

B1 = input( 'B1 = \n' );

%% RF pulse parameters
% to see the slice profile, check out "selective_excitation.m"
% (it is the same pulse)

% nominal flip angle (for B1 == 1 at center of the slice profile)

fa_exc_deg = 90; %input( 'Excitation flip angle [deg] = \n' );
fa_ref_deg = 180; %input( 'Refocusing flip angle [deg] = \n' );

% convert to units, as expected by CoMoTk

fa_exc_rad = fa_exc_deg * pi / 180;
fa_ref_rad = fa_ref_deg * pi / 180;
resolution_um = 1000 * resolution_mm;

%% initialize configuration model (idealized sequence)

cm_cpmg = CoMoTk;
cm_no_cpmg = CoMoTk;

% mandatory tissue parameters

cm_cpmg.R1 = 1 / T1;
cm_cpmg.R2 = 1 / T2;
cm_cpmg.D = D;

cm_no_cpmg.R1 = 1 / T1;
cm_no_cpmg.R2 = 1 / T2;
cm_no_cpmg.D = D;

% further parameters

cm_cpmg.B1 = B1;
cm_no_cpmg.B1 = B1;

% get default options

options = cm_cpmg.options;

options.alloc_n = 1000;  % CoMoTk will allocate more, if needed

options.alloc_d = 2;     % 1 == half echo spacing with crusher
                         % 2 == time between last echo readout and next excitation pulse 
                         %      (includes cpmg spoiler)

options.epsilon = 0;

% set new options

cm_cpmg.options = options;
cm_no_cpmg.options = options;

% start with longitudinal magnetization

cm_cpmg.init_configuration ( [ 0; 0; 1 ] );
cm_no_cpmg.init_configuration ( [ 0; 0; 1 ] );

%% prepare time between end of RF pulse and echo

% unique index

mu_te = 1;

% duration

tau_te = 0.5 * echo_spacing;

% gradient moment == crusher

p_te = [ 2 * pi / resolution_um; 0; 0 ];

%% prepare time between last echo readout and next excitation pulse

% unique index

mu_tr = 2;

% duration

tau_tr = TR - num_echos * echo_spacing;

% gradient moment == Inf (cpmg spoiler, even if D == 0)

p_tr = [ Inf; 0; 0 ];

%% allocate space for results

m_cpmg = zeros( num_echos, num_cycles );
m_no_cpmg = zeros( num_echos, num_cycles );

%% TSE with instantaneous RF pulse

for i = 1 : num_cycles

    % excitation pulse
    
    cm_cpmg.RF( fa_exc_rad, 0 );
    cm_no_cpmg.RF( fa_exc_rad, 0 );

    for j = 1 : num_echos       
        
        % time to refucusing pulse
        
        cm_cpmg.time( mu_te, 'tau', tau_te, 'p', p_te );
        cm_no_cpmg.time( mu_te, 'tau', tau_te, 'p', p_te );
        
        % refocusing pulse
        
        cm_cpmg.RF( fa_ref_rad, 0.5 * pi );
        cm_no_cpmg.RF( fa_ref_rad, 0 );
        
        % time to echo
        
        cm_cpmg.time( mu_te, 'tau', tau_te, 'p', p_te );
        cm_no_cpmg.time( mu_te, 'tau', tau_te, 'p', p_te );
    
        % select the correct configuration:
        % only the zero order configuration along the crusher direction contributes to the voxel signal
        
        b_n = cm_cpmg.find( mu_te, 0 );
        
        % calculate the partial sum
        
        iso = cm_cpmg.isochromat( 0, [], b_n );
        
        % save the echo
        
        m_cpmg( j, i ) = iso.xy;

        % no CPMG
        
        b_n = cm_no_cpmg.find( mu_te, 0 );
        
        % calculate the partial sum
        
        iso = cm_no_cpmg.isochromat( 0, [], b_n );
        
        % save the echo
        
        m_no_cpmg( j, i ) = iso.xy;
        
    end
        
    if ( i < num_cycles )

        % time to next excitation pulse
                
        cm_cpmg.time( mu_tr, 'tau', tau_tr, 'p', p_tr );
        cm_no_cpmg.time( mu_tr, 'tau', tau_tr, 'p', p_tr );
    
    end
        
end

%% Compare the decays of the last cycle

te = ( 1 : num_echos )' .* echo_spacing;

subplot( 1, 2, 1 );
plot( te, abs( m_cpmg( :, end ) ), te, abs( m_no_cpmg( :, end ) ) );
legend( 'CPMG', 'no CPMG' );
title( 'ideal TSE' );

subplot( 1, 2, 2 );
semilogy( te, abs( m_cpmg( :, end ) ), te, abs( m_no_cpmg( :, end ) ) );
legend( 'CPMG', 'no CPMG' );
title( 'ideal TSE' );
