%% GRE sequence
% For GRE sequences (balanced or unbalanced) and a pure tissue,
% this script compares the effect of idealized (instantaneous) versus realistic RF pulses

%% Simulation parameters

par = [];
opt = [];
str = [];

par.ssfp = 'balanced';
par.T1 = 1000;
par.T2 = 100;
par.D = 0;
par.resolution = 1;
par.TR = 5;
par.TE = 2;
par.num_TR = 5;
par.B1 = 1;
par.fa = 50;
par.ph_inc = 180;
par.phd_inc = 0;
par.res = 1;
par.sl_th = 1;
par.t_rf = 1;
par.supp_rf = 100;
par.qual_rf = 3;
par.filt_rf = 'Hamming';
par.epsilon = 1e-4;
par.verbose = 'False';

opt.ssfp = { 'balanced', 'unbalanced' };
opt.T1 = [];
opt.T2 = [];
opt.D = [];
opt.resolution = [];
opt.TR = [];
opt.TE = [];
opt.num_TR = [];
opt.B1 = [];
opt.fa = [];
opt.ph_inc = [];
opt.phd_inc = [];
opt.sl_th = [];
opt.t_rf = [];
opt.supp_rf = [];
opt.qual_rf = [];
opt.filt_rf = { 'Hamming', 'None' };
opt.epsilon = [];
opt.verbose = { 'True', 'False' };

str.ssfp = 'SSFP variant to simulate';
str.T1 = '[ms]';
str.T2 = '[ms]';
str.D = '[um^2/mm] isotropic ADC';
str.resolution = '[mm] in-plane resolution (determines 2*pi crusher, relevant for D ~= 0 only)';
str.TR = '[ms] repetition time';
str.TE = '[ms] echo time (t_rf / 2 <= TE <= TR - t_rf / 2)';
str.num_TR = 'number of TR intervals to simulate';
str.B1 = 'relative B1+';
str.fa = '[deg] flip angle';
str.ph_inc = '[deg] RF phase increment (same as negative off-resonance)' ;
str.phd_inc = '[deg] RF phase difference increment (RF spoiling)';
str.sl_th = '[mm] slice thickness';
str.t_rf = '[ms] RF pulse duration';
str.supp_rf = 'number of RF support points';
str.qual_rf = 'number of SINC pulse zero crossings (on each side)';
str.filt_rf = 'filter RF pulse to reduce wiggles';
str.epsilon = 'discard configuration vectors with L2 norm smaller than this';
str.verbose = 'provide some informal output';

while ( true )
    
    [ par, select_conf ] = sfv( par, opt, str );
    
    if ( select_conf == -1 )
        
        break;
        
    end
    
    %% apply settings
    % calculate phase cycling
    
    ph_deg = zeros( par.num_TR, 1 );
    
    for idx_TR = 2 : par.num_TR
        
        ph_deg( idx_TR ) = ph_deg( idx_TR - 1 ) + par.ph_inc + ( idx_TR - 1 ) * par.phd_inc;
        
    end
    
    % convert to units, as expected by CoMoTk
    
    fa_rad = par.fa * pi / 180;
    ph_rad = ph_deg .* ( pi / 180 );
    sl_th_um = 1000 * par.sl_th;
    
    %% set the RF pulse profile
    % as described in the Handbook of MRI for a SINC pulse +/- Hamming filter
    
    % Number of RF support points
    %
    % as in SLR optimization we approximate the RF pulse as an alternating train of
    %
    % *   par.supp_rf + 1 instantaneous RF pulses, separated by
    % *   par.supp_rf     intervals of constant gradient amplitude
    
    % duration of small intervals
    
    tau_rf = par.t_rf / par.supp_rf;
    
    % support points and flip angles of instantaneous RF pulses (\propto B1 profile)
    
    t = par.qual_rf .* linspace( -1, 1, par.supp_rf + 1 );
    
    if ( isequal( par.filt_rf, 'Hamming' ) )
        
        al = ( 0.54 + 0.46 .* cos( pi .* t ./ par.qual_rf ) ) .* sinc( t );
        
    elseif ( isequal( par.filt_rf, 'None' ) )
        
        al = sinc( t );
        
    end
    
    % normalization: at the center of the slice profile, we need the desired flip angle
    
    al_rad = al .* ( fa_rad / sum( al ) );
    
    % set verbosity
    
    if ( isequal( par.verbose, 'True' ) )
        
        verbose = true;
        
    else
        
        verbose = false;
        
    end
        
    % in-plane crusher gradient moment
    
    p_cru = 2 * pi / ( 1000 * par.resolution );
    
    %% initialize configuration model (idealized sequence)

    cm_ideal = CoMoTk;
    
    % mandatory tissue parameters
    
    cm_ideal.R1 = 1 / par.T1;
    cm_ideal.R2 = 1 / par.T2;
    cm_ideal.D = par.D;
    
    % further parameters
    
    cm_ideal.B1 = par.B1;
    
    % prepare RF pulses
    
    RF_ideal = cell( 1, par.num_TR );
    
    for idx_TR = 1 : par.num_TR

        RF_ideal{ idx_TR }.FlipAngle = fa_rad;
        RF_ideal{ idx_TR }.Phase = ph_rad( idx_TR );
        
    end
    
    % prepare time periods

    RF_to_Echo_ideal = [];
    Echo_to_RF_ideal = [];

    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        RF_to_Echo_ideal.lambda = 1; % unique index
        RF_to_Echo_ideal.tau = par.TE;
        RF_to_Echo_ideal.p = [ 0; 0; 0 ]; % no crusher here ...
        
        Echo_to_RF_ideal.lambda = 2;
        Echo_to_RF_ideal.tau = par.TR - par.TE;
        Echo_to_RF_ideal.p = [ p_cru; 0; 0 ]; % ... but after the echo
        
    elseif ( isequal( par.ssfp, 'balanced' ) )

        RF_to_Echo_ideal.lambda = 1; % unique index
        RF_to_Echo_ideal.tau = par.TR / 2; % for balanced SSFP, we consider a centered echo 
        RF_to_Echo_ideal.p = [ 0; 0; 0 ]; % no crusher
        
        Echo_to_RF_ideal = RF_to_Echo_ideal; % the two intervals are equivalent (duration & (zero) gradient moment)
                
    end
    
    % set options
    
    options = cm_ideal.options;

    options.alloc_n = 1000;  % CoMoTk will allocate more, if needed
    options.alloc_d = 2;     % enough for both variants (unbalanced == 2, balanced == 1)
    options.epsilon = par.epsilon;

    cm_ideal.options = options;
    
    % start with longitudinal magnetization
        
    cm_ideal.init_configuration ( [ 0; 0; 1 ] );
                
    %% initialize configuration model (real sequence)

    cm_real = CoMoTk;
    
    % mandatory tissue parameters
    
    cm_real.R1 = 1 / par.T1;
    cm_real.R2 = 1 / par.T2;
    cm_real.D = par.D;
    
    % further parameters
    
    cm_real.B1 = par.B1;
    
    % approximate real RF pulses as sequence of small instantaneous pulses ...

    RF_real = cell( par.supp_rf + 1, par.num_TR );
    
    for idx_TR = 1 : par.num_TR
        
        for idx_RF = 1 : par.supp_rf + 1
        
            RF_real{ idx_RF, idx_TR }.FlipAngle = al_rad( idx_RF );
            RF_real{ idx_RF, idx_TR }.Phase  = ph_rad( idx_TR );
        
        end
    
    end
        
    % ... separated by small time intervals of constant duration and gradient moment

    DeltaTime = [];
    DeltaTime.lambda = 1;
    DeltaTime.tau = par.t_rf / par.supp_rf;
    
    % the associated gradient moment is defined by RF bandwidth, pulse shape and slice thickness:
    
    % RF pulse bandwidth
    
    f = 2 * par.qual_rf / par.t_rf;
    
    % total moment of slice selection gradient
    
    p_sl = 2 * pi * f * par.t_rf / sl_th_um;
    
    % split into par.supp_rf intervals
    
    DeltaTime.p = [ 0; 0; p_sl / par.supp_rf ];

    % prepare remaining time periods
    
    RF_to_Echo_real = [];
    Echo_to_RF_real = [];
    
    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        RF_to_Echo_real.lambda = 2; % unique index
        RF_to_Echo_real.tau = par.TE - par.t_rf / 2;
        RF_to_Echo_real.p = [ 0; 0; - p_sl / 2 ];
        
        Echo_to_RF_real.lambda = 3;
        Echo_to_RF_real.tau = par.TR - par.TE - par.t_rf / 2;
        Echo_to_RF_real.p = [ p_cru; 0; - p_sl / 2 ];
        
    elseif ( isequal( par.ssfp, 'balanced' ) )
        
        RF_to_Echo_real.lambda = 2; % unique index
        RF_to_Echo_real.tau = ( par.TR - par.t_rf ) / 2; % for balanced SSFP, we consider a centered echo
        RF_to_Echo_real.p = [ 0; 0; - p_sl / 2 ];
        
        Echo_to_RF_real = RF_to_Echo_real; % the two intervals are equivalent (duration & (zero) gradient moment)
        
    end
    
    % set options
    
    options = cm_real.options;
    
    options.alloc_n = 10000;  % CoMoTk will allocate more, if needed
    options.alloc_d = 3;     % enough for both variants (unbalanced == 3, balanced == 2)
    options.epsilon = par.epsilon;
    options.verbose = verbose;
    
    cm_real.options = options;
    
    % start with longitudinal magnetization
    
    cm_real.init_configuration ( [ 0; 0; 1 ] );
                
    %% allocate space for results
        
    m_ideal = zeros( par.num_TR, 1 );
    m_real = zeros( par.num_TR, 1 );

    %% initialize timing
    
    t_id = zeros( par.num_TR, 1 );
    t_re = zeros( par.num_TR, 1 );

    tic;
    t_tmp = toc;

    %%  loop over TR
    
    for idx_TR = 1 : par.num_TR
        
        fprintf( 1, 'i = %d / %d\n', idx_TR, par.num_TR );
                
        %% execute (b)SSFP sequence with ideal, instantaneous pulse ...
        
        if ( idx_TR > 1 ) % not needed before first RF pulse
            
            cm_ideal.time( Echo_to_RF_ideal );
            
        end
        
        cm_ideal.RF( RF_ideal{ idx_TR } );
        
        cm_ideal.time( RF_to_Echo_ideal );
          
        %% ... and extract the result
        
        % first we select the configurations
        % this depends on the situation

        if ( isequal( par.ssfp, 'unbalanced' ) )
                        
            % we assume a conventional FID sequence with crusher AFTER the echo
            % i.e. in the interval 'Echo_to_RF'

            % the associated configuration must be zero
            % since the signal from nonzero orders is (hopefully) dephased 
            
            select_conf = [];
            select_conf.b_n = cm_ideal.find( Echo_to_RF_ideal.lambda, 0 );
            
        elseif ( isequal( par.ssfp, 'balanced' ) )
            
            % balanced SSFP has no crusher gradients
            % in absence of gradient dephasing, all occupied configurations contribute

            select_conf = [];
            select_conf.b_n = cm_ideal.b_n;
            
        end

        % calculate the partial sum
        
        res = cm_ideal.sum( select_conf );
        
        % save the echo
        
        m_ideal( idx_TR ) = res.xy;
            
        % time spent for the ideal part
        
        t_id( idx_TR ) = - t_tmp;
        t_tmp = toc;
        t_id( idx_TR ) = t_id( idx_TR ) + t_tmp;        
        
        %% execute (b)SSFP with real RF pulse ...
        
        if ( idx_TR > 1 ) % not needed before first RF pulse
            
            cm_real.time( Echo_to_RF_real );
            
        end

        for idx_RF = 1 : par.supp_rf
        
            cm_real.RF( RF_real{ idx_RF, idx_TR } );

            cm_real.time( DeltaTime );
        
        end
        
        cm_real.RF( RF_real{ par.supp_rf + 1, idx_TR } );

        cm_real.time( RF_to_Echo_real );

        %% ... and extract the result
            
        if ( isequal( par.ssfp, 'unbalanced' ) )
                        
            % we assume a conventional FID sequence with crusher AFTER the echo
            % i.e. in the interval 'Echo_to_RF'

            % like in the real case, the associated configuration must be zero
            % since the signal from nonzero orders is (hopefully) dephased 
            
            select_conf = [];
            select_conf.b_n = cm_real.find( Echo_to_RF_real.lambda, 0 );
                                    
        elseif ( isequal( par.ssfp, 'balanced' ) )
            
            % balanced SSFP has no crusher gradients
            % in absence of gradient dephasing, all occupied configurations contribute

            select_conf = [];
            select_conf.b_n = cm_real.b_n;

        end
        
        % in both cases and unlike the real case, we *further* need to integrate over the slice profile
            
        select_conf.b_n = select_conf.b_n & ...
            reshape( abs( cm_real.p_n( 3, : ) ) < 0.5 * abs( DeltaTime.p( 3 ) ), size( select_conf.b_n ) );
        
                % calculate the partial sum
        
        res = cm_real.sum( select_conf );
        
        % save the echo
        
        m_real( idx_TR ) = res.xy;
            
        % time spent for the real part
        
        t_re( idx_TR ) = - t_tmp;
        t_tmp = toc;
        t_re( idx_TR ) = t_re( idx_TR ) + t_tmp;        
        
        fprintf( 1, 'ideal RF = %9.3f sec\n', t_id( idx_TR ) );
        fprintf( 1, 'real RF  = %9.3f sec\n', t_re( idx_TR ) );
        
    end
    
    fprintf( 1, 'ideal RF =\t %9.3f sec total\n', sum( t_id ) );
    fprintf( 1, 'real RF  =\t %9.3f sec total\n', sum( t_re ) );
    
    %% Show results

    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        te = par.TE + ( 0 : par.num_TR - 1 )' .* par.TR;
        
    else
        
        te = par.TR / 2 + ( 0 : par.num_TR - 1 )' .* par.TR;
        
    end
   
    sc = abs( m_ideal( 1 ) / m_real( 1 ) );

    subplot( 1, 1, 1 );
    plot( te, abs( m_ideal ), te, sc .* abs( m_real ) );
    legend( 'ideal', 'real' );

    if ( isequal( par.ssfp, 'unbalanced' ) )
        
        title( 'SSFP' );
        
    elseif ( isequal( par.ssfp, 'balanced' ) )
        
        title( 'bSSFP' );
        
    end
    
end
