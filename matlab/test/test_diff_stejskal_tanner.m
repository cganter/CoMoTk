%% Calculates diffusion damping for Stejskal-Tanner gradients
% for a sequence as described in the documentation:
%
% * SE sequence with instantaneous 90° and 180° pulses
% * delta : duration of rectangular gradients
% * Delta : relative displacement of gradients
% 
% the user specifies the duration delta <= tau together with the 
% starting times 0 <= t_1, t_2 <= tau - delta
% this defines Delta = tau + t2 - t1
%
% for b = 0 : 0.25 : 2, suitable gradient moments are calculated with
%
% p = sqrt( b / ( Delta - delta / 3 ) )
% 
% the gradient shape of interval j = 1,2 is defined by the two parameters
%
% s_j = 2 * p * ( tau - t_j - delta / 2 ) / tau
%
% S_j = p^2 * ( tau - t_j - 2 * delta / 3 ) / tau

%% Simulation parameters
% (by default, all set by the user)

% Tissue properties

T1 = input( 'T1 [ms] = \n' );
T2 = input( 'T2 [ms] = \n' );
D = input( 'D [um^2/ms] =\n' ); % water approximately corresponds to a value of 3

% Sequence and measurement settings

TE = input( 'TE [ms] = \n' );
tau = TE / 2;

% specify duration of gradient pulse

delta = -1;

str = [ 'choose a delta between 0 and ', num2str( tau ), ':\n' ];

while( delta < 0 || delta > tau )

    delta = input( str );
    
end

% starting times

t_s = - ones( 2, 1 );

% specify first starting time

str = [ 'choose first starting time between 0 and ', num2str( tau - delta ), ':\n' ];

while( t_s( 1 ) < 0 || t_s( 1 ) > tau - delta )

    t_s( 1 ) = input( str );
    
end

% specify second starting time

str = [ 'choose second starting time between 0 and ', num2str( tau - delta ), ':\n' ];

while( t_s( 2 ) < 0 || t_s( 2 ) > tau - delta )

    t_s( 2 ) = input( str );
    
end

% this defines Delta

Delta = tau + t_s( 2 ) - t_s( 1 );

% the b-values of interest

b = 0 : 0.25 : 2;

% the corresponding gradient moments

p_ = sqrt( b ./ ( Delta - delta ./ 3 ) );

% and their shapes

s_ = 2 .* p_ .* ( tau - t_s - 0.5 .* delta ) ./ tau;

S_ = p_.^2 .* ( tau - t_s - 2 .* delta ./ 3 ) ./ tau;

% allocate space for results

m_xy = zeros( size( b ) );

%% perform the simulation for each b-value

for i = 1 : length( b )

    %% initialize configuration model (idealized sequence)

    cm = CoMoTk;   % constant gradient assumed implicitly

    % mandatory tissue parameters

    cm.R1 = 1 / T1;
    cm.R2 = 1 / T2;
    cm.D = D;

    % get default options

    options = cm.options;

    options.alloc_n = 10;  % more than needed...
    options.alloc_d = 1;   % only one time interval
    options.epsilon = 0;   % do not discard any modes

    % set new options

    cm.options = options;

    % start with longitudinal magnetization

    cm.init_configuration ( [ 0; 0; 1 ] );

    %% prepare time between two RF pulses

    % unique index

    mu_time = 1;

    % gradient moment

    p = [ p_( i ); 0; 0 ];

    % excitation pulse

    cm.RF( pi / 2, 0 );

    % first time period
    % gradient shape

    s = [ s_( 1, i ); 0; 0; S_( 1, i ) ];

    cm.time( mu_time, 'tau', tau, 'p', p, 's', s );

    % refocusing pulse with CPMG condition (not that it would matter here...)

    cm.RF( pi, pi / 2 );

    % second time period
    % gradient shape

    s = [ s_( 2, i ); 0; 0; S_( 2, i ) ];

    cm.time( mu_time, 'tau', tau, 'p', p, 's', s );

    % read echo signal

    m_xy( i ) = sqrt( 2 ) * cm.m( 1, cm.find( mu_time, 0 ) );

end

%% show results

subplot( 1, 2, 1 );
plot( b, abs( m_xy ), b, exp( - b * D - TE / T2 ) );
legend( 'simulation', 'theory' );

subplot( 1, 2, 2 );
semilogy( b, abs( m_xy ), b, exp( - b * D - TE / T2 ) );
legend( 'simulation', 'theory' );
