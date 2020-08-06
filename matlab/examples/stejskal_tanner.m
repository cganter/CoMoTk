%% Calculates diffusion damping for Stejskal-Tanner gradients
% for a sequence as described in the manuscript:
%
% * SE sequence with instantaneous 90° and 180° pulses
% * delta : duration of rectangular gradients
% * Delta : relative displacement of gradients
% 
% the user specifies the duration delta <= tau together with the 
% starting times 0 <= t1, t2 <= tau - delta
% this defines Delta = tau + t2 - t1
%
% for an user-defined equidistant set of b-values suitable gradient moments are calculated with
%
% p = sqrt( b / ( Delta - delta / 3 ) )
% 
% the gradient shape of interval j = 1,2 is defined by the two parameters
%
% s_j = 2 * p * ( tau - tj - delta / 2 ) / tau
% S_j = p^2 * ( tau - tj - 2 * delta / 3 ) / tau
%
% it is shown that the such defined sequence, after evaluation with the configuration model,
% generates the expected diffusion damping proportional to exp(- b * D )

par = [];
opt = [];
str = [];

par.T1 = 100;
par.T2 = 10;
par.D = 3;
par.tau = 1;
par.delta = 0.5;
par.t1 = 0.25;
par.t2 = 0.25;
par.min_b = 0;
par.max_b = 2;
par.num_b = 10;

opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.tau = [];
opt.delta = [];
opt.t1 = [];
opt.t2 = [];
opt.min_b = [];
opt.max_b = [];
opt.num_b = [];

str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/ms] isotropic ADC';
str.tau = '[ms] (== TE/2)';
str.delta = '[ms] gradient duration (0 < delta < tau';
str.t1 = '[ms] start of first gradient (0 <= t1 <= tau - delta)';
str.t2 = '[ms] start of second gradient (0 <= t2 <= tau - delta)';
str.min_b = '[mm/um^2] minimal b-value'; 
str.max_b = '[mm/um^2] maximal b-value'; 
str.num_b = 'number of b-value to be simulated'; 

while ( true )
    
    [ par, sel ] = set_field_values( par, opt, str );
    
    if ( sel == -1 )
        
        break;
        
    end
    
    %% Simulation parameters
    % (by default, all set by the user)
        
    % Tissue properties
    
    T1 = par.T1;
    T2 = par.T2;
    D = par.D;
    
    % Sequence and measurement settings
    
    tau = par.tau;
    TE = 2 * tau;
    
    % Diffusion gradient duration and placement

    delta = par.delta;
    t_s = [ par.t1; par.t2 ];
    Delta = tau + t_s( 2 ) - t_s( 1 );

    if ( sel == 0 && ( ...
            delta <= 0 || ...
            min( t_s ) < 0 || ...
            max( t_s ) + delta > tau ) )
        
        fprintf( 1, 'Inconsistent gradient timing. Please correct.\n' );
       
        continue;
        
    end
    
    % b-values
    
    b = linspace( par.min_b, par.max_b, par.num_b );
    
    % the corresponding gradient moments
    
    p_ = sqrt( b ./ ( Delta - delta ./ 3 ) );
    
    % and their shapes (according to the integrals in the manuscript)
    
    s_ = 2 .* p_ .* ( tau - t_s - 0.5 .* delta ) ./ tau;
    S_ = p_.^2 .* ( tau - t_s - 2 .* delta ./ 3 ) ./ tau;
    
    % allocate space for results
    
    m_xy = zeros( size( b ) );
    
    %% prepare parameters
    % RF pulses 
    RF_exc.FlipAngle = pi / 2;
    RF_exc.Phase = 0;

    RF_refoc.FlipAngle = pi;
    RF_refoc.Phase = pi / 2;
  
    % time intervals (constant properties)
    T_1.lambda = 1;
    T_1.tau = tau;
    
    T_2.lambda = 1;
    T_2.tau = tau;
    
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
        options.alloc_d = 1;   % only one time interval and instantaneous RF pulses
        options.epsilon = 0;   % do not discard any modes (== optimal accuracy)
        
        % set new options
        
        cm.options = options;
        
        %% comlete parameters for time intervals (gradient moment and shape)
        T_1.p = [ p_( i ); 0; 0 ];   
        T_1.s = [ s_( 1, i ); 0; 0 ];
        T_1.S = S_( 1, i );

        T_2.p = [ p_( i ); 0; 0 ];
        T_2.s = [ s_( 2, i ); 0; 0 ];
        T_2.S = S_( 2, i );

        %% run sequence
        cm.RF( RF_exc );
        cm.time( T_1 );
        cm.RF( RF_refoc );
        cm.time( T_2 );
        
        %% read echo signal
        % select zero order configuration
        conf_0 = [];
        conf_0.b_n = cm.find( T_1.lambda, 0 );
            
        % calculate the partial sum
        res = cm.sum( conf_0 );
        
        % save the echo
        m_xy( i ) = res.xy;
                
    end
    
    %% show results
    
    subplot( 1, 2, 1 );
    plot( b, abs( m_xy ), 'o', b, exp( - b * D - TE / T2 ) );
    legend( 'simulation', 'theory' );
    
    subplot( 1, 2, 2 );
    semilogy( b, abs( m_xy ), 'o', b, exp( - b * D - TE / T2 ) );
    legend( 'simulation', 'theory' );
    
end
