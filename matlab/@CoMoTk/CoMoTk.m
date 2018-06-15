classdef CoMoTk < matlab.mixin.Copyable
    % CoMoTk  Configuration Model Toolkit
    %
    % Features:
    %
    % - optional: diffusion effects
    % - optional: derivatives with respect to R1, R2, B1 and D
    %
    %
    % Carl Ganter 2018
    
    properties ( Constant )
        
        % frequently used
        
        sqrt_2 = sqrt( 2 );
        sqrt_2_inv = sqrt( 0.5 );
        
    end
    
    properties ( Dependent )
        
        % global options
        
        options;
        
        % ALL dependent quantities (except for the optional s)
        % have to be set explicitly prior to simulation. Otherwise an error is returned.
        
        %% tissue properties
        % relaxation rates
        
        R1;                                                                % [ 1 / ms ]
        R2;                                                                % [ 1 / ms ]
        
        % diffusion sensitivity (scalar or tensor)
        
        D;                                                                 % [ um^2 / ms ]  i.e. D( Water ) is approximately equal to 3
        
        %% multiple tissues
        % relative weight
        
        w;
        
        % delta omega
        
        dom;                                                             % [ radian / ms ]
        
        % R2^\prime
        
        R2p;                                                             % [ 1 / ms ]
        
        %% measurement conditions
        % relative B1
        
        B1;
        
    end
    
    properties ( SetAccess = private )
        
        %% global options (to be set prior to calling CoMoTk.init_configuration)
        % the purpose of the alloc_* variables is to increase the numerical performance by minimizing
        % -> recalculation of identical quantities (alloc_RF)
        % -> reallocation of memory (alloc_d, alloc_n)
        % the set values define
        % -> the initial state
        % -> and the minimial reallocation increment, if a limit is reached
        
        options_priv = struct( ...
            'alloc_d', 1, ...                                            % dimension of configuration model
            'alloc_n', 1000, ...                                         % number of occupied states
            'epsilon', 0, ...                                         % discard configuration vectors with a squared norm smaller than this
            'rapid_meltdown', true, ...                            % faster, but more memory needed
            'verbose', false, ...                                      % various status messages
            'debug', false ...                                         % store state for debugging purposes
            );
        
        % actual values
        
        alloc_d = [];
        alloc_n = [];
        epsilon = [];
        rapid_meltdown = [];
        verbose = [];
        debug = [];
        
        %% information about the actual state
        
        m = [];                                                            % configuration vector
        
        b_mu = [];                                                         % occupied dimensions
        b_n = [];                                                          % occupied configurations (total)
        
        mu = [];                                                           % indices of dimensions
        n = [];                                                            % configuration orders in each dimension
        
        n_tissues = 0;		                                               % number of tissues
        d = 0;                                                             % == sum( cm.b_mu ) : dimension of configuration model
        n_conf = 0;                                                        % == sum( cm.b_n ) : occupied configurations (including n == 0)
        null_idx = 1;                                                      % index, corresponding to n == 0 (remains untouched)
        
        %% mandatory tissue properties
        % cf. dependent counterparts
        
        R1_priv = [];
        R2_priv = [];
        D_priv = [];
        
        %% optional tissue properties and measurement conditions
        % cf. dependent counterparts
        
        w_priv = 1;
        dom_priv = 0;
        R2p_priv = 0;
        B1_priv = 1;
        
        %% sequence properties
        % time per interval
        
        b_tau = [];
        tau = [];                                                          % [ ms ]
        
        % gradient moment (for each time period)
        % unit: [ 1 / um ] (spatial frequency)
        % definition: gamma * integral( dt * g( t ) )
        % with gamma == 2 * pi * 42.576 / ( mT * ms )
        % and g( t ) the gradient in [ mT / um ]
        
        b_p = [];
        p = [];                                                            % [ 1 / um ]
        
        % integrals for arbitrary (but otherwise constant) gradient shapes
        % see documentation for the precise definition
        
        b_s = [];
        s = [];                                                            % [ 1 / um ]   for s( 1 : 3, : )
        % [ 1 / um^2 ] for s( 4, : )
        
        %% time and gradient moment, associated with configuration order
        % updated after every time interval
        
        b_tau_n = [];
        tau_n = [];
        
        b_p_n = [];
        p_n = [];
        
        %% derivatives
        % handles of derivatives to calculate
        
        dR1 = [];
        dR2 = [];
        dD = [];
        dB1 = [];
        dFlipAngle = [];
        dPhase = [];
        dtau = [];
        dp = [];
        ds = [];
        
        % associated numbers ...
        
        len_dR1 = 0;
        len_dR2 = 0;
        len_dD = 0;
        len_dB1 = 0;
        len_dFlipAngle = 0;
        len_dPhase = 0;
        len_dtau = 0;
        len_dp = 0;
        len_ds = 0;
        
        % ... and derivatives of configuration vector
        
        dm_dR1 = [];
        dm_dR2 = [];
        dm_dD = [];
        dm_dB1 = [];
        dm_dFlipAngle = [];
        dm_dPhase = [];
        dm_dtau = [];
        dm_dp = [];
        dm_ds = [];
        
        %% auxiliary variables to speed up calculations
        % bookkeeping of occupied states
        
        b_up_free = [];
        b_do_free = [];
        
        idx_up = [];
        idx_do = [];
        
        idx_up_new = [];
        idx_do_new = [];
        
        % pure relaxation
        
        b_E = [];
        E = [];
        
        % scalar products in diffusion exponent and their derivatives
        
        b_pz = [];
        pz = [];
        dpz_dp = [];
        
        b_pxy = [];
        pxy = [];
        dpxy_dp = [];
        dpxy_ds = [];
        
        % relaxation with diffusion and derivatives
        
        b_EFn = [];
        EFn = [];
        
        b_dEFn_dD = [];
        dEFn_dD = [];
        
        b_dEFn_dtau = [];
        dEFn_dtau = [];
        
        b_dEFn_dp = [];
        dEFn_dp = [];
        
        b_dEFn_ds = [];
        dEFn_ds = [];
        
        % debugging
        
        log = [];
        
    end
    
    methods
        
        %% constructor
        
        function cm = CoMoTk ( )
            % CONSTRUCTOR
            % no arguments, no action needed
            
        end
        
        %% set / get methods (mandatory variables)
        
        function set.R1 ( cm, R1 )
            
            if ( ~isempty( cm.m ) )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.init_configuration.' );
                
            end
            
            if ( min( R1 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R1 );
                cm.w = ones( 1, cm.n_tissues );
                
            elseif ( length( R1 ) ~= cm.n_tissues )
                
                error( 'length( R1 ) must be equal to number of tissues.' );
                
            end
            
            cm.R1_priv = R1( : )';                                         % row vector
            
        end
        
        function res = get.R1 ( cm )
            
            res = cm.R1_priv;
            
        end
        
        function set.R2 ( cm, R2 )
            
            if ( ~isempty( cm.m ) )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.init_configuration.' );
                
            end
            
            if ( min( R2 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R2 );
                cm.w = ones( 1, cm.n_tissues );
                
            elseif ( length( R2 ) ~= cm.n_tissues )
                
                error( 'length( R2 ) must be equal to number of tissues.' );
                
            end
            
            cm.R2_priv = R2( : )';                                         % row vector
            
        end
        
        function res = get.R2 ( cm )
            
            res = cm.R2_priv;
            
        end
        
        function set.D ( cm, D )
            
            if ( ~isempty( cm.m ) )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.init_configuration.' );
                
            end
            
            if ( min( D ) < 0 )
                
                error( 'Diffusion constant must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( D );
                cm.w = ones( 1, cm.n_tissues );
                
            elseif ( length( D ) ~= cm.n_tissues )
                
                error( 'length( D ) must be equal to number of tissues.' );
                
            end
            
            cm.D_priv = D( : )';                                           % row vector
            
        end
        
        function res = get.D ( cm )
            
            res = cm.D_priv;
            
        end
        
        %% set / get methods (optional variables)
        
        % relative B1 (default == 1)
        
        function set.B1 ( cm, B1 )
            
            if ( ~isempty( cm.m ) )
                
                error( 'B1 can only be set prior to calling CoMoTk.init_configuration.' );
                
            end
            
            if ( B1 < 0 )
                
                error( 'B1 must not be negative.' );
                
            end
            
            cm.B1_priv = B1;
            
        end
        
        function res = get.B1 ( cm )
            
            res = cm.B1_priv;
            
        end
        
        % individual weights (== proton density) of complex tissues
        % (default == 1)
        
        function set.w ( cm, w )
            
            if ( ~isempty( cm.m ) )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.init_configuration.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( w );
                
            elseif ( length( w ) ~= cm.n_tissues )
                
                error( 'length( w ) must be equal to number of tissues.' );
                
            end
            
            cm.w_priv = reshape( w, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.w ( cm )
            
            res = cm.w_priv;
            
        end
        
        % relative frequencies
        % (default == 0)
        
        function set.dom ( cm, dom )
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( dom );
                
            elseif ( length( dom ) ~= cm.n_tissues )
                
                error( 'length( dom ) must be equal to number of tissues.' );
                
            end
            
            cm.dom_priv = reshape( dom, [ 1, 1, cm.n_tissues ] );
            
        end
        
        function res = get.dom ( cm )
            
            res = cm.dom_priv;
            
        end
        
        % R2^\prime
        % (default == 0)
        
        function set.R2p ( cm, R2p )
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R2p );
                
            elseif ( length( R2p ) ~= cm.n_tissues )
                
                error( 'length( R2p ) must be equal to number of tissues.' );
                
            end
            
            cm.R2p_priv = reshape( R2p, [ 1, 1, cm.n_tissues ] );
            
        end
        
        function res = get.R2p ( cm )
            
            res = cm.R2p_priv;
            
        end
        
        % options
        
        function set.options ( cm, options )
            
            if ( ~isempty( cm.m ) )
                
                error( 'Options can only be set prior to calling CoMoTk.init_configuration.' );
                
            end
            
            cm.options_priv = options;
            
        end
        
        function res = get.options ( cm )
            
            res = cm.options_priv;
            
        end
        
        %% prepare magnetization
        
        function init_configuration ( cm, m )
            % INIT_CONFIGURATION  Initializes configuration
            %
            % IN
            %
            % m : Initial real magnetization vector ( row or column ) with
            %     m( 1 ) == m_x
            %     m( 2 ) == m_y
            %     m( 3 ) == m_z
            %
            % Assumption: only zero order configuration occupied
            
            cm.check_mandatory_parameters;
            
            if ( isempty( m ) )  % initialize with proton densities
                
                m = zeros( 3, cm.n_tissues );
                m( 3, : ) = cm.w;
                
            elseif ( size( m, 1 ) ~= 3 || size( m, 2 ) ~= cm.n_tissues )
                
                error( 'Magnetization has wrong dimension.' );
                
            end
            
            % take initial settings from options
            
            cm.alloc_d = cm.options.alloc_d;
            cm.alloc_n = cm.options.alloc_n;
            cm.epsilon = cm.options.epsilon;
            cm.rapid_meltdown = cm.options.rapid_meltdown;
            cm.verbose = cm.options.verbose;
            cm.debug = cm.options.debug;
            
            % generate initial configuration vector
            
            cm.m = zeros( 3, cm.n_tissues, cm.alloc_n );
            
            cm.m( 1, :, cm.null_idx ) = ( m( 1, : ) + 1i * m( 2, : ) ) * CoMoTk.sqrt_2_inv;
            cm.m( 2, :, cm.null_idx ) = m( 3, : );
            cm.m( 3, :, cm.null_idx ) = conj( cm.m( 1, :, cm.null_idx ) );
            
            % allocate space
            
            cm.b_mu = false( 1, cm.alloc_d );
            cm.b_n = false( cm.alloc_n, 1 );
            
            % cm.mu = zeros( 1, cm.alloc_d, 'uint32' );
            % cm.n = zeros( cm.alloc_n, cm.alloc_d, 'int32' );
            
            cm.mu = zeros( 1, cm.alloc_d );
            cm.n = zeros( cm.alloc_n, cm.alloc_d );
            
            cm.b_up_free = true( cm.alloc_n, cm.alloc_d );
            cm.b_do_free = true( cm.alloc_n, cm.alloc_d );
            
            %             cm.idx_up = zeros( cm.alloc_n, cm.alloc_d, 'uint32' );
            %             cm.idx_do = zeros( cm.alloc_n, cm.alloc_d, 'uint32' );
            %
            %             cm.idx_up_new = zeros( cm.alloc_n, 1, 'uint32' );
            %             cm.idx_do_new = zeros( cm.alloc_n, 1, 'uint32' );
            
            cm.idx_up = zeros( cm.alloc_n, cm.alloc_d );
            cm.idx_do = zeros( cm.alloc_n, cm.alloc_d );
            
            cm.idx_up_new = zeros( cm.alloc_n, 1 );
            cm.idx_do_new = zeros( cm.alloc_n, 1 );
            
            % only the zeroth configuration is occupied
            
            cm.b_n( cm.null_idx ) = true;
            cm.n_conf = 1;
            
            % initialize times and gradients
            
            cm.b_tau = false( cm.alloc_d, 1 );
            cm.tau = zeros( cm.alloc_d, 1 );
            
            cm.b_p = false( cm.alloc_d, 1 );
            cm.p = zeros( 3, cm.alloc_d );
            
            cm.b_s = false( cm.alloc_d, 1 );
            cm.s = zeros( 4, cm.alloc_d );
            
            cm.b_tau_n = false( cm.alloc_n, 1 );
            cm.tau_n = zeros( cm.alloc_n, 1 );

            cm.b_p_n = false( cm.alloc_n, 1 );
            cm.p_n = zeros( 3, cm.alloc_n );
            
            % pure relaxation
            
            cm.b_E = false( cm.alloc_d, 1 );
            cm.E = zeros( 2, cm.n_tissues, cm.alloc_d );
            
        end
        
        %% define, which derivatives should be calculated
        
        function set_derivatives ( cm, varargin )
            % SET_DERIVATIVES  Define the set of derivatives, to be calculated
            %
            % IN
            %
            % varargin : pairs (X,Y), where
            %            X is from the set { 'R1', 'R2', 'D', 'B1', 'FlipAngle', 'Phase', 'tau', 'p', 's' }
            %            and Y is an integer
            %            exception: for X == 'B1', the argument Y is not given
            
            if ( isempty( cm.m ) )
                
                error( 'Derivatives must be set after calling CoMoTk.init_configuration.' );
                
            end
            
            nVarargs = length( varargin );
            
            i = 1;
            
            while( i <= nVarargs )
                
                if ( isequal( varargin{ i }, 'R1' ) )
                    
                    cm.dR1 = varargin{ i + 1 };
                    cm.len_dR1 = length( cm.dR1 );
                    cm.dm_dR1 = zeros( 3, cm.len_dR1, cm.alloc_n );
                    
                elseif ( isequal( varargin{ i }, 'R2' ) )
                    
                    cm.dR2 = varargin{ i + 1 };
                    cm.len_dR2 = length( cm.dR2 );
                    cm.dm_dR2 = zeros( 3, cm.len_dR2, cm.alloc_n );
                    
                elseif ( isequal( varargin{ i }, 'D' ) )
                    
                    cm.dD = varargin{ i + 1 };
                    cm.len_dD = length( cm.dD );
                    cm.dm_dD =   zeros( 3, cm.len_dD, cm.alloc_n );
                    cm.dEFn_dD = zeros( 3, cm.len_dD, cm.alloc_n, cm.alloc_d );
                    
                elseif ( isequal( varargin{ i }, 'B1' ) )
                    
                    cm.dB1 = 1;
                    cm.len_dB1 = 1;
                    cm.dm_dB1 = zeros( 3, cm.n_tissues, cm.alloc_n );
                    i = i - 1;
                    
                elseif ( isequal( varargin{ i }, 'FlipAngle' ) )
                    
                    cm.dFlipAngle = varargin{ i + 1 };
                    cm.len_dFlipAngle = length( cm.dFlipAngle );
                    cm.dm_dFlipAngle = zeros( 3, cm.n_tissues, cm.alloc_n, cm.len_dFlipAngle );
                    
                elseif ( isequal( varargin{ i }, 'Phase' ) )
                    
                    cm.dPhase = varargin{ i + 1 };
                    cm.len_dPhase = length( cm.dPhase );
                    cm.dm_dPhase = zeros( 3, cm.n_tissues, cm.alloc_n, cm.len_dPhase );
                    
                elseif ( isequal( varargin{ i }, 'tau' ) )
                    
                    cm.dtau = varargin{ i + 1 };
                    cm.len_dtau = length( cm.dtau );
                    cm.dm_dtau   = zeros( 3, cm.n_tissues, cm.alloc_n, cm.len_dtau );
                    cm.dEFn_dtau = zeros( 3, cm.n_tissues, cm.alloc_n, cm.len_dtau );
                    
                elseif ( isequal( varargin{ i }, 'p' ) )
                    
                    cm.dp = varargin{ i + 1 };
                    cm.len_dp = length( cm.dp );
                    cm.dm_dp   = zeros( 3, cm.n_tissues, cm.alloc_n, 3, cm.len_dp );
                    cm.dpz_dp  = zeros( 1,            1, cm.alloc_n, 3, cm.len_dp );
                    cm.dpxy_dp = zeros( 2,            1, cm.alloc_n, 3, cm.len_dp, cm.alloc_d );
                    cm.dEFn_dp = zeros( 3, cm.n_tissues, cm.alloc_n, 3, cm.len_dp, cm.alloc_d );
                    
                elseif ( isequal( varargin{ i }, 's' ) )
                    
                    cm.ds = varargin{ i + 1 };
                    cm.len_ds = length( cm.ds );
                    cm.dm_ds   = zeros( 3, cm.n_tissues, cm.alloc_n, 4, cm.len_ds );
                    cm.dpxy_ds = zeros( 2,            1, cm.alloc_n, 4, cm.len_ds, cm.alloc_d );
                    cm.dEFn_ds = zeros( 2, cm.n_tissues, cm.alloc_n, 4, cm.len_ds, cm.alloc_d );
                    
                else
                    
                    error( 'Invalid argument.' );
                    
                end
                
                i = i + 2;
                
            end
            
        end
        
        %% spin dynamics
        
        function RF ( cm, FlipAngle, Phase, varargin )
            % RF  Applies an instantaneous RF pulse
            %
            % IN
            %
            % alpha    : flip angle [rad] (scalar)
            % phase    : Phase of RF pulse [rad] (scalar)
            % varargin : optional ( tag, handle ) pairs
            %
            % format of varargin:
            %
            % tag = 'Pulse' | 'FlipAngle' | 'Phase'
            % handle = positive integer
            %
            % meaning of arguments:
            %
            % 'Pulse' : If the same RF pulse is executed more than once, it makes sense
            %           to store the rotation matrix (including possible derivatives)
            %           at the first execution for reuse in later calls.
            %           Each pulse is uniquely defined by its handle.
            %
            % OUT
            %
            % out   : state immediately after instantaneous rotation
            % block : information, which can be used to speed up next call
            %         via method CoMoTk.run (of course only, if identical
            %         RF pulse is repeated)
            %
            % All phases simply refer to some given fixed coordinate
            % system and are not defined relative to the previous RF pulse.
            
            % check, whether all variables have been set.
            
            if ( isempty( cm.m ) )
                
                error( 'Incomplete specification of configuration model.' );
                
            end
            
            % this is our actual flip angle
            
            ActualFlipAngle = cm.B1 * FlipAngle;
            
            % no action required, if no rotation is performed
            
            if ( mod( ActualFlipAngle, 2 * pi ) == 0 )
                
                return;
                
            end
            
            % check for derivative handles
            
            handle_FlipAngle = 0;
            handle_Phase = 0;
            
            nVarargs = length( varargin );
            
            for i = 1 : 2 : nVarargs
                
                if ( isequal( varargin{ i }, 'FlipAngle' ) && ...
                        i < nVarargs && ...
                        isnumeric( varargin{ i + 1 } ) )
                    
                    handle_FlipAngle = varargin{ i + 1 };
                    
                elseif ( isequal( varargin{ i }, 'Phase' ) && ...
                        i < nVarargs && ...
                        isnumeric( varargin{ i + 1 } ) )
                    
                    handle_Phase = varargin{ i + 1 };
                    
                else
                    
                    error( 'Invalid handle argument.' );
                    
                end
                
            end
            
            % calculate the rotation matrix
            
            c_al = cos( ActualFlipAngle );
            s_al = sin( ActualFlipAngle );
            
            ei_ph = cos( Phase ) + 1i * sin( Phase );
            
            RotMat = zeros( 3 );
            
            RotMat( 1, 1 ) = 0.5 * ( 1 + c_al );
            RotMat( 1, 2 ) = - 1i * s_al * ei_ph * CoMoTk.sqrt_2_inv;
            RotMat( 1, 3 ) = 0.5 * ( 1 - c_al ) * ei_ph^2;
            RotMat( 2, 1 ) = - conj( RotMat( 1, 2 ) );
            RotMat( 2, 2 ) = c_al;
            RotMat( 2, 3 ) = - RotMat( 1, 2 );
            RotMat( 3, 1 ) = conj( RotMat( 1, 3 ) );
            RotMat( 3, 2 ) = - RotMat( 2, 1 );
            RotMat( 3, 3 ) = RotMat( 1, 1 );
            
            %% calculate derivatives, if necessary
            % This has to be done (at least for B1, FlipAngle and Phase) prior
            % to updating the configuration vector cm.m itself.
            % The rotation matrix does not depend on R1, R2, D, tau, p and s
            % which simplifes the associated recursions.
            
            % R1 derivative(s)
            
            if ( cm.len_dR1 > 0 )
                
                cm.dm_dR1( :, :, cm.b_n ) = reshape( ...
                    RotMat * reshape( cm.dm_dR1( :, :, cm.b_n ), [ 3, cm.len_dR1 * cm.n_conf ] ), ...
                    [ 3, cm.len_dR1, cm.n_conf ] );
                
            end
            
            % R2 derivative(s)
            
            if ( cm.len_dR2 > 0 )
                
                cm.dm_dR2( :, :, cm.b_n ) = reshape( ...
                    RotMat * reshape( cm.dm_dR2( :, :, cm.b_n ), [ 3, cm.len_dR2 * cm.n_conf ] ), ...
                    [ 3, cm.len_dR2, cm.n_conf ] );
                
            end
            
            % D derivative(s)
            
            if ( cm.len_dD > 0 )
                
                cm.dm_dD( :, :, cm.b_n ) = reshape( ...
                    RotMat * reshape( cm.dm_dD( :, :, cm.b_n ), [ 3, cm.len_dD * cm.n_conf ] ), ...
                    [ 3, cm.len_dD, cm.n_conf ] );
                
            end
            
            % B1 derivative
            
            if ( cm.len_dB1 > 0 )
                
                dc_al = - s_al * FlipAngle;
                ds_al = c_al * FlipAngle;
                
                % calculate the rotation matrix
                
                dRotMat_dB1 = zeros( 3 );
                
                dRotMat_dB1( 1, 1 ) = 0.5 * dc_al;
                dRotMat_dB1( 1, 2 ) = - 1i * ds_al * ei_ph * CoMoTk.sqrt_2_inv;
                dRotMat_dB1( 1, 3 ) = - 0.5 * dc_al * ei_ph^2;
                dRotMat_dB1( 2, 1 ) = - conj( dRotMat_dB1( 1, 2 ) );
                dRotMat_dB1( 2, 2 ) = dc_al;
                dRotMat_dB1( 2, 3 ) = - dRotMat_dB1( 1, 2 );
                dRotMat_dB1( 3, 1 ) = conj( dRotMat_dB1( 1, 3 ) );
                dRotMat_dB1( 3, 2 ) = - dRotMat_dB1( 2, 1 );
                dRotMat_dB1( 3, 3 ) = dRotMat_dB1( 1, 1 );
                
                cm.dm_dB1( :, :, cm.b_n ) = reshape( ...
                    RotMat * reshape( cm.dm_dB1( :, :, cm.b_n ), [ 3, cm.n_tissues * cm.n_conf ] ) + ...
                    dRotMat_dB1 * reshape( cm.m( :, :, cm.b_n ), [ 3, cm.n_tissues * cm.n_conf ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf ] );
                
            end
            
            % FlipAngle derivative(s)
            
            if ( cm.len_dFlipAngle > 0 )
                
                cm.dm_dFlipAngle( :, :, cm.b_n, : ) = reshape( ...
                    RotMat * reshape( cm.dm_dFlipAngle( :, :, cm.b_n, : ), [ 3, cm.n_tissues * cm.n_conf * cm.len_dFlipAngle ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf, cm.len_dFlipAngle ] );
                
                [ ~, i ] = find( cm.dFlipAngle == handle_FlipAngle );
                
                if ( ~isempty( i ) )
                    
                    dc_al = - s_al * cm.B1;
                    ds_al = c_al * cm.B1;
                    
                    % calculate the rotation matrix
                    
                    dRotMat_dFlipAngle = zeros( 3 );
                    
                    dRotMat_dFlipAngle( 1, 1 ) = 0.5 * dc_al;
                    dRotMat_dFlipAngle( 1, 2 ) = - 1i * ds_al * ei_ph * CoMoTk.sqrt_2_inv;
                    dRotMat_dFlipAngle( 1, 3 ) = - 0.5 * dc_al * ei_ph^2;
                    dRotMat_dFlipAngle( 2, 1 ) = - conj( dRotMat_dFlipAngle( 1, 2 ) );
                    dRotMat_dFlipAngle( 2, 2 ) = dc_al;
                    dRotMat_dFlipAngle( 2, 3 ) = - dRotMat_dFlipAngle( 1, 2 );
                    dRotMat_dFlipAngle( 3, 1 ) = conj( dRotMat_dFlipAngle( 1, 3 ) );
                    dRotMat_dFlipAngle( 3, 2 ) = - dRotMat_dFlipAngle( 2, 1 );
                    dRotMat_dFlipAngle( 3, 3 ) = dRotMat_dFlipAngle( 1, 1 );
                    
                    cm.dm_dFlipAngle( :, :, cm.b_n, i ) = cm.dm_dFlipAngle( :, :, cm.b_n, i ) + ...
                        reshape( ...
                        dRotMat_dFlipAngle * reshape( cm.m( :, :, cm.b_n ), [ 3, cm.n_tissues * cm.n_conf ] ), ...
                        [ 3, cm.n_tissues, cm.n_conf ] );
                    
                end
                
            end
            
            % Phase derivative(s)
            
            if ( cm.len_dPhase > 0 )
                
                cm.dm_dPhase( :, :, cm.b_n, : ) = reshape( ...
                    RotMat * reshape( cm.dm_dPhase( :, :, cm.b_n, : ), [ 3, cm.n_tissues * cm.n_conf * cm.len_dPhase ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf, cm.len_dPhase ] );
                
                [ ~, i ] = find( cm.dPhase == handle_Phase );
                
                if ( ~isempty( i ) )
                    
                    % calculate the rotation matrix
                    
                    dRotMat_dPhase = zeros( 3 );
                    
                    dRotMat_dPhase( 1, 2 ) = 1i * Rotmat( 1, 2 );
                    dRotMat_dPhase( 1, 3 ) = 2 * 1i * Rotmat( 1, 3 );
                    dRotMat_dPhase( 2, 1 ) = - conj( dRotMat_dPhase( 1, 2 ) );
                    dRotMat_dPhase( 2, 3 ) = - dRotMat_dPhase( 1, 2 );
                    dRotMat_dPhase( 3, 1 ) = conj( dRotMat_dPhase( 1, 3 ) );
                    dRotMat_dPhase( 3, 2 ) = - dRotMat_dPhase( 2, 1 );
                    
                    cm.dm_dPhase( :, :, cm.b_n, i ) = cm.dm_dPhase( :, :, cm.b_n, i ) + ...
                        reshape( ...
                        dRotMat_dPhase * reshape( cm.m( :, :, cm.b_n ), [ 3, cm.n_tissues * cm.n_conf ] ), ...
                        [ 3, cm.n_tissues, cm.n_conf ] );
                    
                end
                
            end
            
            % tau derivative(s)
            
            if( cm.len_dtau > 0 )
                
                cm.dm_dtau( :, :, cm.b_n, : ) = reshape( ...
                    RotMat * reshape( cm.dm_dtau( :, :, cm.b_n, : ), [ 3, cm.n_tissues * cm.n_conf * cm.len_dtau ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf, cm.len_dtau ] );
                
            end
            
            % p derivative(s)
            
            if( cm.len_dp > 0 )
                
                cm.dm_dp( :, :, cm.b_n, :, : ) = reshape( ...
                    RotMat * reshape( cm.dm_dp( :, :, cm.b_n, :, : ), [ 3, cm.n_tissues * cm.n_conf * 3 * cm.len_dp ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf, 3, cm.len_dp ] );
                
            end
            
            % s derivative(s)
            
            if( cm.len_ds > 0 )
                
                cm.dm_ds( :, :, cm.b_n, :, : ) = reshape( ...
                    RotMat * reshape( cm.dm_ds( :, :, cm.b_n, :, : ), [ 3, cm.n_tissues * cm.n_conf * 4 * cm.len_ds ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf, 4, cm.len_ds ] );
                
            end
            
            %% finally we update the configuration cm.m
            
            % the reshapes are required by the matrix multiplication
            
            cm.m( :, :, cm.b_n ) = reshape( ...
                RotMat * reshape( cm.m( :, :, cm.b_n ), [ 3, cm.n_tissues * cm.n_conf ] ), ...
                [ 3, cm.n_tissues, cm.n_conf ] );
            
        end
        
        function b_is_ok = time ( cm, handle_time, varargin )
            % TIME  time evolution during period T( handle_time ) (without RF pulse)
            %
            % IN
            %
            % in    : prior state (configuration vector and optional derivatives)
            %
            % OUT
            %
            % out   : state immediately after time period
            % block : information, which can be used to speed up next call
            %         via method CoMoTk.run
            %
            % Different from the instantaneous pulse, the states of adjacent
            % configuration orders mix during this period.
            
            %% update dimensions and indices
            
            b_is_ok = true;
            
            if ( cm.options.debug && ~cm.check_state( 'before pre_time' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            [ idx_time, b_n_new, b_up_free_time, b_do_free_time ] = cm.pre_time( handle_time );
            
            if ( cm.options.debug && ~cm.check_state( 'after pre_time' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            %% check optional arguments ( tau, p, s )
            
            b_new_mu = ~cm.b_tau( idx_time );
            
            nVarargs = length( varargin );
            
            i = 1;
            
            while( i <= nVarargs )
                
                if ( isequal( varargin{ i }, 'p' ) )
                    
                    if ( b_new_mu )
                        
                        cm.b_p( idx_time ) = true;
                        cm.p( :, idx_time ) = varargin{ i + 1 }( : );
                        
                    else
                        
                        if ( max( abs( cm.p( :, idx_time ) - varargin{ i + 1 }( : ) ) ) > 0 )

                            error( 'Gradient moment p must not change.' );
                        
                        end
                            
                    end
                    
                elseif ( isequal( varargin{ i }, 's' ) )
                    
                    % Derivatives with respect to s make only sense for those intervals,
                    % for which the shape is not variable.
                    % This is checked by the following statement.
                    
                    if ( b_new_mu )
                        
                        cm.b_s( idx_time ) = true;
                        
                    else
                        
                        [ ~, i_ ] = find( cm.ds == handle_time );
                        
                        if ( ~isempty( i_ ) && ...
                              max( abs( cm.s( :, i_ ) - varargin{ i + 1 }( : ) ) ) > 0 )
                            
                            error( 'Derivative of variable shape makes no sense.' );
                            
                        end
                        
                    end
                    
                    cm.s( :, idx_time ) = varargin{ i + 1 }( : );
                    
                elseif ( isequal( varargin{ i }, 'tau' ) )
                    
                    if ( b_new_mu )
                        
                        cm.b_tau( idx_time ) = true;
                        cm.tau( idx_time ) = varargin{ i + 1 };
                        
                        if ( cm.tau( idx_time ) <= 0 )
                            
                            error( 'Time interval must be positive.' );
                            
                        end
                        
                    elseif ( cm.tau( idx_time ) ~= varargin{ i + 1 } )
                        
                        error( 'Time interval tau must not change.' );
                        
                    end
                    
                else
                    
                    error( 'Invalid argument.' );
                    
                end
                
                i = i + 2;
                
            end
            
            % at this point, a nonzero (positive) time should have been set
            
            if ( ~cm.b_tau( idx_time ) )
                
                error( 'Failed to specify time interval.' );
                
            end
            
            % rule out a few invalid combinations
            
            if( cm.b_s( idx_time ) )
                
                if ( ~cm.b_p( idx_time ) )
                    
                    error( 'Definition of a gradient shape without specifying its moment is not allowed.' );
                    
                end
                
                if( max( abs( cm.p( :, idx_time ) ) ) == Inf && cm.s( 4, idx_time ) ~= Inf )
                    
                    error( 'Inconsistent gradient specification.' );
                    
                end
                
            end
            
            %% update relaxation matrices
            
            cm.update_relaxation( idx_time );
            
            %% calculate recursive update of derivatives, if necessary
            % At least for R1, R2, D, tau, p and s, this has to be done prior
            % to updating the configuration vector cm.m itself.
            % The recursion coefficients do not depend on B1, FlipAngle and Phase
            % which simplifes the associated recursions.
            
            % R1 derivative(s)
            
            if ( cm.len_dR1 > 0 )
                
                if ( isempty( cm.EFn ) )
                    
                    cm.dm_dR1( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                        cm.E( 2, cm.dR1, idx_time ) .* cm.dm_dR1( 1, :, cm.b_n );
                    
                    cm.dm_dR1( 2, :, cm.b_n ) = cm.E( 1, cm.dR1, idx_time ) .* ...
                        ( cm.dm_dR1( 2, :, cm.b_n ) - cm.tau( idx_time ) .* cm.m( 2, cm.dR1, cm.b_n ) );
                    
                    cm.dm_dR1( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                        cm.E( 2, cm.dR1, idx_time ) .* cm.dm_dR1( 3, :, cm.b_n );
                    
                else
                    
                    cm.dm_dR1( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                        cm.EFn( 2, cm.dR1, cm.b_n, idx_time ) .* cm.dm_dR1( 1, :, cm.b_n );
                    
                    cm.dm_dR1( 2, :, cm.b_n ) = cm.EFn( 2, :, cm.b_n, idx_time ) .* ...
                        ( cm.dm_dR1( 2, :, cm.b_n ) - cm.tau( idx_time ) .* cm.m( 2, cm.dR1, cm.b_n ) ) ;
                    
                    cm.dm_dR1( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                        cm.EFn( 3, cm.dR1, cm.b_n, idx_time ) .* cm.dm_dR1( 3, :, cm.b_n );
                    
                end
                
                cm.dm_dR1( 2, :, cm.null_idx ) = cm.dm_dR1( 2, :, cm.null_idx ) + ...
                    cm.w( cm.dR1 ) .* cm.tau( idx_time ) .* cm.E( 1, cm.dR1, idx_time );
                
            end
            
            % R2 derivative(s)
            
            if ( cm.len_dR2 > 0 )
                
                if ( isempty( cm.EFn ) )
                    
                    cm.dm_dR2( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = cm.E( 2, cm.dR2, idx_time ) .* ...
                        ( cm.dm_dR2( 1, :, cm.b_n ) - cm.tau( idx_time ) .* cm.m( 1, cm.dR2, cm.b_n ) );
                    
                    cm.dm_dR2( 2, :, cm.b_n ) = ...
                        cm.E( 1, cm.dR2, idx_time ) .* cm.dm_dR2( 2, :, cm.b_n );
                    
                    cm.dm_dR2( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = cm.E( 2, cm.dR2, idx_time ) .* ...
                        ( cm.dm_dR2( 3, :, cm.b_n ) - cm.tau( idx_time ) .* cm.m( 3, cm.dR2, cm.b_n ) );
                    
                else
                    
                    cm.dm_dR2( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = cm.EFn( 1, cm.dR2, cm.b_n, idx_time ) .* ...
                        ( cm.dm_dR2( 1, :, cm.b_n ) - cm.tau( idx_time ) .* cm.m( 1, cm.dR2, cm.b_n ) );
                    
                    cm.dm_dR2( 2, :, cm.b_n ) = ...
                        cm.EFn( 2, cm.dR2, cm.b_n, idx_time ) .* cm.dm_dR2( 2, :, cm.b_n );
                    
                    cm.dm_dR2( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = cm.EFn( 3, cm.dR2, cm.b_n, idx_time ) .* ...
                        ( cm.dm_dR2( 3, :, cm.b_n ) - cm.tau( idx_time ) .* cm.m( 3, cm.dR2, cm.b_n ) );
                    
                end
                
            end
            
            % D derivative(s)
            
            if ( cm.len_dD > 0 )
                
                cm.dm_dD( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                    cm.EFn( 1, cm.dD, cm.b_n, idx_time ) .* cm.dm_dD( 1, :, cm.b_n ) + ...
                    cm.dEFn_dD( 1, cm.dD, cm.b_n, idx_time ) .* cm.m( 1, cm.dD, cm.b_n );
                
                cm.dm_dD( 2, :, cm.b_n ) = ...
                    cm.EFn( 2, cm.dD, cm.b_n, idx_time ) .* cm.dm_dD( 2, :, cm.b_n ) + ...
                    cm.dEFn_dD( 2, cm.dD, cm.b_n, idx_time ) .* cm.m( 2, cm.dD, cm.b_n );
                
                cm.dm_dD( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                    cm.EFn( 3, cm.dD, cm.b_n, idx_time ) .* cm.dm_dD( 3, :, cm.b_n ) + ...
                    cm.dEFn_dD( 3, cm.dD, cm.b_n, idx_time ) .* cm.m( 3, cm.dD, cm.b_n );
                
            end
            
            % B1 derivative
            
            if ( cm.len_dB1 > 0 )
                
                if ( isempty( cm.EFn ) )
                    
                    cm.dm_dB1( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dB1( 1, :, cm.b_n );
                    
                    cm.dm_dB1( 2, :, cm.b_n ) = ...
                        cm.E( 1, :, idx_time ) .* cm.dm_dB1( 2, :, cm.b_n );
                    
                    cm.dm_dB1( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dB1( 3, :, cm.b_n );
                    
                else
                    
                    cm.dm_dB1( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                        cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.dm_dB1( 1, :, cm.b_n );
                    
                    cm.dm_dB1( 2, :, cm.b_n ) = ...
                        cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.dm_dB1( 2, :, cm.b_n );
                    
                    cm.dm_dB1( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                        cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.dm_dB1( 3, :, cm.b_n );
                    
                end
                
            end
            
            % FlipAngle derivative(s)
            
            if ( cm.len_dFlipAngle > 0 )
                
                if ( isempty( cm.EFn ) )
                    
                    cm.dm_dFlipAngle( 1, :, cm.idx_up( cm.b_n, idx_time ), : ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dFlipAngle( 1, :, cm.b_n, : );
                    
                    cm.dm_dFlipAngle( 2, :, cm.b_n, : ) = ...
                        cm.E( 1, :, idx_time ) .* cm.dm_dFlipAngle( 2, :, cm.b_n, : );
                    
                    cm.dm_dFlipAngle( 3, :, cm.idx_do( cm.b_n, idx_time ), : ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dFlipAngle( 3, :, cm.b_n, : );
                    
                else
                    
                    cm.dm_dFlipAngle( 1, :, cm.idx_up( cm.b_n, idx_time ), : ) = ...
                        cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.dm_dFlipAngle( 1, :, cm.b_n, : );
                    
                    cm.dm_dFlipAngle( 2, :, cm.b_n, : ) = ...
                        cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.dm_dFlipAngle( 2, :, cm.b_n, : );
                    
                    cm.dm_dFlipAngle( 3, :, cm.idx_do( cm.b_n, idx_time ), : ) = ...
                        cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.dm_dFlipAngle( 3, :, cm.b_n, : );
                    
                end
                
            end
            
            % Phase derivative(s)
            
            if ( cm.len_dPhase > 0 )
                
                if ( isempty( cm.EFn ) )
                    
                    cm.dm_dPhase( 1, :, cm.idx_up( cm.b_n, idx_time ), : ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dPhase( 1, :, cm.b_n, : );
                    
                    cm.dm_dPhase( 2, :, cm.b_n, : ) = ...
                        cm.E( 1, :, idx_time ) .* cm.dm_dPhase( 2, :, cm.b_n, : );
                    
                    cm.dm_dPhase( 3, :, cm.idx_do( cm.b_n, idx_time ), : ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dPhase( 3, :, cm.b_n, : );
                    
                else
                    
                    cm.dm_dPhase( 1, :, cm.idx_up( cm.b_n, idx_time ), : ) = ...
                        cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.dm_dPhase( 1, :, cm.b_n, : );
                    
                    cm.dm_dPhase( 2, :, cm.b_n, : ) = ...
                        cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.dm_dPhase( 2, :, cm.b_n, : );
                    
                    cm.dm_dPhase( 3, :, cm.idx_do( cm.b_n, idx_time ), : ) = ...
                        cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.dm_dPhase( 3, :, cm.b_n, : );
                    
                end
                
            end
            
            % tau derivative(s)
            
            if ( cm.len_dtau > 0 )
                
                if ( isempty( cm.EFn ) )
                    
                    cm.dm_dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), : ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dtau( 1, :, cm.b_n, : );
                    
                    cm.dm_dtau( 2, :, cm.b_n, : ) = ...
                        cm.E( 1, :, idx_time ) .* cm.dm_dtau( 2, :, cm.b_n, : );
                    
                    cm.dm_dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), : ) = ...
                        cm.E( 2, :, idx_time ) .* cm.dm_dtau( 3, :, cm.b_n, : );
                    
                else
                    
                    cm.dm_dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), : ) = ...
                        cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.dm_dtau( 1, :, cm.b_n, : );
                    
                    cm.dm_dtau( 2, :, cm.b_n, : ) = ...
                        cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.dm_dtau( 2, :, cm.b_n, : );
                    
                    cm.dm_dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), : ) = ...
                        cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.dm_dtau( 3, :, cm.b_n, : );
                    
                end
                
                [ ~, i ] = find( cm.dtau == handle_time );
                
                if ( ~isempty( i ) )
                    
                    if ( isempty( cm.EFn ) )
                        
                        cm.dm_dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), i ) = cm.dm_dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), i ) - ...
                            cm.R2 .* cm.E( 2, :, idx_time ) .* cm.m( 1, :, cm.b_n );
                        
                        cm.dm_dtau( 2, :, cm.b_n, i ) = cm.dm_dtau( 2, :, cm.b_n, i ) - ...
                            cm.R1 .* cm.E( 1, :, idx_time ) .* cm.m( 2, :, cm.b_n );
                        
                        cm.dm_dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), i ) = cm.dm_dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), i ) - ...
                            cm.R2 .* cm.E( 2, :, idx_time ) .* cm.m( 3, :, cm.b_n );
                        
                    else
                        
                        cm.dm_dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), i ) = cm.dm_dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), i ) + ...
                            cm.dEFn_dtau( 1, :, cm.b_n, idx_time ) .* cm.m( 1, :, cm.b_n );
                        
                        cm.dm_dtau( 2, :, cm.b_n, i ) = cm.dm_dtau( 2, :, cm.b_n, i ) + ...
                            cm.dEFn_dtau( 2, :, cm.b_n, idx_time ) .* cm.m( 2, :, cm.b_n );
                        
                        cm.dm_dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), i ) = cm.dm_dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), i ) + ...
                            cm.dEFn_dtau( 3, :, cm.b_n, idx_time ) .* cm.m( 3, :, cm.b_n );
                        
                    end
                    
                    cm.dm_dtau( 2, :, cm.null_idx, i ) = cm.dm_dtau( 2, :, cm.null_idx, i ) + cm.w .* cm.R1 .* cm.E( 1, :, idx_time );
                    
                end
                
            end
            
            % p derivative(s)
            
            if ( cm.len_dp > 0 )
                
                cm.dm_dp( 1, :, cm.idx_up( cm.b_n, idx_time ), :, : ) = ...
                    cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.dm_dp( 1, :, cm.b_n, :, : ) + ...
                    cm.dEFn_dp( 1, :, cm.b_n, :, :, idx_time ) .* cm.m( 1, :, cm.b_n );
                
                cm.dm_dp( 2, :, cm.b_n, :, : ) = ...
                    cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.dm_dp( 2, :, cm.b_n, :, : ) + ...
                    cm.dEFn_dp( 2, :, cm.b_n, :, :, idx_time ) .* cm.m( 2, :, cm.b_n );
                
                cm.dm_dp( 3, :, cm.idx_do( cm.b_n, idx_time ), :, : ) = ...
                    cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.dm_dp( 3, :, cm.b_n, :, : ) + ...
                    cm.dEFn_dp( 3, :, cm.b_n, :, :, idx_time ) .* cm.m( 3, :, cm.b_n );
                
            end
            
            % s derivative(s)
            
            if ( cm.len_ds > 0 )
                
                cm.dm_ds( 1, :, cm.idx_up( cm.b_n, idx_time ), :, : ) = ...
                    cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.dm_ds( 1, :, cm.b_n, :, : );
                
                cm.dm_ds( 2, :, cm.b_n, :, : ) = ...
                    cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.dm_ds( 2, :, cm.b_n, :, : );
                
                cm.dm_ds( 3, :, cm.idx_do( cm.b_n, idx_time ), :, : ) = ...
                    cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.dm_ds( 3, :, cm.b_n, :, : );
                
                [ ~, i ] = find( cm.ds == handle_time );
                
                if ( ~isempty( i ) )
                    
                    cm.dm_ds( 1, :, cm.idx_up( cm.b_n, idx_time ), :, i ) = cm.dm_ds( 1, :, cm.idx_up( cm.b_n, idx_time ), :, i ) + ...
                        cm.dEFn_ds( 1, :, cm.b_n, :, i, idx_time ) .* cm.m( 1, :, cm.b_n );
                    
                    cm.dm_ds( 3, :, cm.idx_do( cm.b_n, idx_time ), :, i ) = cm.dm_ds( 3, :, cm.idx_do( cm.b_n, idx_time ), :, i ) + ...
                        cm.dEFn_ds( 2, :, cm.b_n, :, i, idx_time ) .* cm.m( 3, :, cm.b_n );
                    
                end
                
            end
            
            %% finally we update the configuration cm.m
            
            % relaxation
            
            if ( isempty( cm.EFn ) )
                
                cm.m( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                    cm.E( 2, :, idx_time ) .* cm.m( 1, :, cm.b_n );
                
                cm.m( 2, :, cm.b_n ) = ...
                    cm.E( 1, :, idx_time ) .* cm.m( 2, :, cm.b_n );
                
                cm.m( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                    cm.E( 2, :, idx_time ) .* cm.m( 3, :, cm.b_n );
                
                cm.m( 1, :, b_do_free_time ) = 0;
                cm.m( 3, :, b_up_free_time ) = 0;
                
            else
                
                cm.m( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                    cm.EFn( 1, :, cm.b_n, idx_time ) .* cm.m( 1, :, cm.b_n );
                
                cm.m( 2, :, cm.b_n ) = ...
                    cm.EFn( 2, :, cm.b_n, idx_time ) .* cm.m( 2, :, cm.b_n );
                
                cm.m( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                    cm.EFn( 3, :, cm.b_n, idx_time ) .* cm.m( 3, :, cm.b_n );
                
                cm.m( 1, :, b_do_free_time ) = 0;
                cm.m( 3, :, b_up_free_time ) = 0;
                
            end
            
            % repolarization (here, the relative proton densities come into play)
            
            cm.m( 2, :, cm.null_idx ) = cm.m( 2, :, cm.null_idx ) + cm.w .* ( 1 - cm.E( 1, :, idx_time ) );
            
            %% perform final bookkeeping
            
            cm.post_time( b_n_new );
            
            if ( cm.options.debug && ~cm.check_state( 'after post_time' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            %% remove negligible configurations
            
            if( ~cm.meltdown )
                
                b_is_ok = false;
                return;
                
            end
            
            if ( cm.options.debug && ~cm.check_state( 'after meltdown' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            %% update configuration dependent time and gradient
            
            cm.update_tau_n;
            cm.update_p_n;
            
        end
        
        %% get results
        
        function b_n = find ( cm, mu, n ) % <<find>>
            
            if ( length( mu ) ~= length( n ) )
                
                error( 'mu and n must have same length.' );
                
            end
            
            % check, which of the supplied mu are actually alive (i.e. stored)
            
            b_mu_2d = cm.b_mu & ( mu( : ) == cm.mu );
            
            % non-existing mu must not have non-zero configuration orders
            
            b_mu_found = any( b_mu_2d, 2 );
            
            if ( any( n( ~b_mu_found ) ) )
                
                b_n = [];
                return;
                
            end
            
            n_mu_found = sum( b_mu_found );
            
            % for the rest, determine the indices in storage
            
            [ idx_mu, ~ ] = find( b_mu_2d' );
            
            % we can now determine b_n
            
            b_n = cm.b_n & any( cm.n( :, idx_mu ) == reshape( n( b_mu_found ), [ 1, n_mu_found ] ), 2 );
            
        end
        
        function iso = isochromat ( cm, om, x, b_n )
            
            if ( isempty( b_n ) )
                
                b_n = cm.b_n;
                
            end
            
            n_c = sum( b_n );
            
            % off-resonance and susceptibility effects
            
            osg = ( 1i .* ( om + cm.dom ) - cm.R2p ) .* cm.tau( cm.b_mu ).';
            
            % plus gradients (optional)
            
            if ( ~isempty( x ) )
                
                osg = osg - 1i .* sum( cm.p( :, cm.b_mu ) .* x( : ) );
                
            end
            
            % weighted by configuration order and added together
            
            n_osg = sum( cm.n( b_n, cm.b_mu ) .* osg, 2 );
            
            % some reordering and reshaping
            
            if ( isempty( n_osg ) )
                
                n_osg = 0;
                
            else
                
                n_osg = reshape( reshape( n_osg, [ n_c, cm.n_tissues ] ).', [ 1, cm.n_tissues, n_c ] );
                
            end
            
            % phase factor, scaled by susceptibility effects
            
            e_n_osg = exp( n_osg );
            
            % calculate the isochromat
            
            iso.xy = CoMoTk.sqrt_2 .* sum( sum( e_n_osg .* cm.m( 1, :, b_n ), 3 ) );
            iso.z = sum( sum ( e_n_osg .* cm.m( 2, :, b_n ), 3 ) );  % should be real...
            
            % and the calculated derivatives
            
            if ( cm.len_dR1 > 0 )
                
                iso.dm_dR1.xy = CoMoTk.sqrt_2 .* sum( e_n_osg( 1, cm.dR1, : ) .* cm.dm_dR1( 1, :, b_n ), 3 );
                iso.dm_dR1.z = sum( e_n_osg( 1, cm.dR1, : ) .* cm.dm_dR1( 2, :, b_n ), 3 );
                
                iso.dm_dR1.xy = reshape( iso.dm_dR1.xy, [ 1, cm.len_dR1 ] );
                iso.dm_dR1.z = reshape( iso.dm_dR1.z, [ 1, cm.len_dR1 ] );
                
            end
            
            if ( cm.len_dR2 > 0 )
                
                iso.dm_dR2.xy = CoMoTk.sqrt_2 .* sum( e_n_osg( 1, cm.dR2, : ) .* cm.dm_dR2( 1, :, b_n ), 3 );
                iso.dm_dR2.z = sum( e_n_osg( 1, cm.dR2, : ) .* cm.dm_dR2( 2, :, b_n ), 3 );
                
                iso.dm_dR2.xy = reshape( iso.dm_dR2.xy, [ 1, cm.len_dR2 ] );
                iso.dm_dR2.z = reshape( iso.dm_dR2.z, [ 1, cm.len_dR2 ] );
                
            end
            
            if ( cm.len_dD > 0 )
                
                iso.dm_dD.xy = CoMoTk.sqrt_2 .* sum( e_n_osg( 1, cm.dD, : ) .* cm.dm_dD( 1, :, b_n ), 3 );
                iso.dm_dD.z = sum( e_n_osg( 1, cm.dD, : ) .* cm.dm_dD( 2, :, b_n ), 3 );
                
                iso.dm_dD.xy = reshape( iso.dm_dD.xy, [ 1, cm.len_dD ] );
                iso.dm_dD.z = reshape( iso.dm_dD.z, [ 1, cm.len_dD ] );
                
            end
            
            if ( cm.len_dB1 > 0 )
                
                iso.dm_dB1.xy = CoMoTk.sqrt_2 .* sum( sum( e_n_osg .* cm.dm_dB1( 1, :, b_n ), 3 ) );
                iso.dm_dB1.z = sum( sum( e_n_osg .* cm.dm_dB1( 2, :, b_n ), 3 ) );
                
            end
            
            if ( cm.len_dFlipAngle > 0 )
                
                iso.dm_dFlipAngle.xy = CoMoTk.sqrt_2 .* ...
                    sum( sum( e_n_osg .* cm.dm_dFlipAngle( 1, :, b_n, : ), 3 ), 2 );
                iso.dm_dFlipAngle.z = sum( sum ( e_n_osg .* cm.dm_dFlipAngle( 2, :, b_n, : ), 3 ), 2 );
                
                iso.dm_dFlipAngle.xy = reshape( iso.dm_dFlipAngle.xy, [ 1, cm.len_dFlipAngle ] );
                iso.dm_dFlipAngle.z = reshape( iso.dm_dFlipAngle.z, [ 1, cm.len_dFlipAngle ] );
                
            end
            
            if ( cm.len_dPhase > 0 )
                
                iso.dm_dPhase.xy = CoMoTk.sqrt_2 .* sum( sum( e_n_osg .* cm.dm_dPhase( 1, :, b_n, : ), 3 ), 2 );
                iso.dm_dPhase.z = sum( sum ( e_n_osg .* cm.dm_dPhase( 2, :, b_n, : ), 3 ), 2 );
                
                iso.dm_dPhase.xy = reshape( iso.dm_dPhase.xy, [ 1, cm.len_dPhase ] );
                iso.dm_dPhase.z = reshape( iso.dm_dPhase.z, [ 1, cm.len_dPhase ] );
                
            end
            
            if ( cm.len_dp > 0 )
                
                iso.dm_dp.xy = CoMoTk.sqrt_2 .* sum( sum( e_n_osg .* cm.dm_dp( 1, :, b_n, :, : ), 3 ), 2 );
                iso.dm_dp.z = sum( sum ( e_n_osg .* cm.dm_dp( 2, :, b_n, :, : ), 3 ), 2 );
                
                iso.dm_dp.xy = reshape( iso.dm_dp.xy, [ 3, cm.len_dp ] );
                iso.dm_dp.z = reshape( iso.dm_dp.z, [ 3, cm.len_dp ] );
                
            end
            
            if ( cm.len_ds > 0 )
                
                iso.dm_ds.xy = CoMoTk.sqrt_2 .* sum( sum( e_n_osg .* cm.dm_ds( 1, :, b_n, :, : ), 3 ), 2 );
                iso.dm_ds.z = sum( sum ( e_n_osg .* cm.dm_ds( 2, :, b_n, :, : ), 3 ), 2 );
                
                iso.dm_ds.xy = reshape( iso.dm_ds.xy, [ 4, cm.len_ds ] );
                iso.dm_ds.z = reshape( iso.dm_ds.z, [ 4, cm.len_ds ] );
                
            end
            
        end
        
        %% various checks for debugging purposes
        
        function b_is_ok = check_state ( cm, context )
            
            fprintf( 1, 'check_state: %s\n', context );
            
            b_is_ok = true;
            
            if ( cm.options.debug )
                
                if ( isempty( cm.log )  )
                    
                    cm.log.count = 1;
                    
                else
                    
                    cm.log.count = cm.log.count + 1;
                    
                end
                
                cm.log.context{ cm.log.count } = context;
                
                cm.log.d{ cm.log.count } = cm.d;
                cm.log.n_conf{ cm.log.count } = cm.n_conf;
                
                cm.log.b_n{ cm.log.count } = cm.b_n;
                cm.log.n{ cm.log.count } = cm.n;
                
                cm.log.b_mu{ cm.log.count } = cm.b_mu;
                cm.log.mu{ cm.log.count } = cm.mu;
                
                cm.log.b_up_free{ cm.log.count } = cm.b_up_free;
                cm.log.b_do_free{ cm.log.count } = cm.b_do_free;
                
                cm.log.idx_up{ cm.log.count } = cm.idx_up;
                cm.log.idx_do{ cm.log.count } = cm.idx_do;
                
            end
            
            if ( cm.options.verbose )
                
                b_up_occ = cm.b_n & ~cm.b_up_free & cm.b_mu;
                b_do_occ = cm.b_n & ~cm.b_do_free & cm.b_mu;
                
                [ idx_up_occ, c_up ] = find( b_up_occ );
                [ idx_do_occ, c_do ] = find( b_do_occ );
                
                if ( ( ~isempty( cm.idx_up( b_up_occ ) ) && min( cm.idx_up( b_up_occ ) ) == 0 ) || ...
                        ( ~isempty( cm.idx_do( b_do_occ ) ) && min( cm.idx_do( b_do_occ ) ) == 0 ) )
                    
                    b_is_ok = false;
                    return;
                    
                end
                
                idx_up_diff = cm.idx_do( sub2ind( [ cm.alloc_n, cm.alloc_d ], cm.idx_up( b_up_occ ), c_up ) ) - idx_up_occ;
                b_ok = idx_up_diff == 0;
                
                if ( cm.options.debug )
                    
                    cm.log.idx_up_diff{ cm.log.count } = idx_up_diff;
                    cm.log.b_idx_up_diff{ cm.log.count } = b_ok;
                    
                end
                
%                fprintf( 1, 'idx_up_diff: %d / %d ok.\n', sum( b_ok ), length( idx_up_diff ) );
                
                if ( sum( b_ok ) ~= length( idx_up_diff ) )
                    
                    b_is_ok = false;
                    
                end
                
                idx_do_diff = cm.idx_up( sub2ind( [ cm.alloc_n, cm.alloc_d ], cm.idx_do( b_do_occ ), c_do ) ) - idx_do_occ;
                b_ok = idx_do_diff == 0;
                
                if ( cm.options.debug )
                    
                    cm.log.idx_do_diff{ cm.log.count } = idx_do_diff;
                    cm.log.b_idx_do_diff{ cm.log.count } = b_ok;
                    
                end
                
%                fprintf( 1, 'idx_do_diff: %d / %d ok.\n', sum( b_ok ), length( idx_do_diff ) );
                
                if ( sum( b_ok ) ~= length( idx_do_diff ) )
                    
                    b_is_ok = false;
                    
                end
                
                n_diff = cm.n( sub2ind( [ cm.alloc_n, cm.alloc_d ], cm.idx_up( b_up_occ ), c_up ) ) - cm.n( b_up_occ );
                b_ok = n_diff == 1;
                
                if ( cm.options.debug )
                    
                    cm.log.n_diff{ cm.log.count } = n_diff;
                    cm.log.b_n_diff{ cm.log.count } = b_ok;
                    
                end
                
%                fprintf( 1, 'n_diff: %d / %d ok.\n', sum( b_ok ), length( n_diff ) );
                
                if ( sum( b_ok ) ~= length( n_diff ) )
                    
                    b_is_ok = false;
                    
                end
                
            end
            
        end
        
    end
    
    methods ( Access = private )
        
        function update_tau_n ( cm )
            
            b_todo = cm.b_n & ~cm.b_tau_n;
            
            n_todo = sum( b_todo );
            
            if ( n_todo > 0 )
                
                cm.tau_n( b_todo ) = ...
                    sum( cm.n( b_todo, cm.b_mu ) .* reshape( cm.tau( cm.b_mu ), [ 1, cm.d ] ), 2 );
                
                cm.b_tau_n = cm.b_n;
                
            end
            
        end
        
        function update_p_n ( cm )
            
            b_todo = cm.b_n & ~cm.b_p_n;
            
            n_todo = sum( b_todo );
            
            if ( n_todo > 0 )
                
                tmp = reshape( cm.n( b_todo, cm.b_mu ), [ 1, n_todo, cm.d ] ) .* ...
                    reshape( cm.p( :, cm.b_mu ), [ 3, 1, cm.d ] );
                
                % correct for 0 * Inf == NaN
                
                tmp( isnan( tmp ) ) = 0;
                
                cm.p_n( :, b_todo ) = sum( tmp, 3 );
                
                cm.b_p_n = cm.b_n;
                
            end
            
        end
        
        function check_mandatory_parameters ( cm )
            
            if ( isempty( cm.R1 ) || ...
                    isempty( cm.R2 ) || ...
                    isempty( cm.D ) )
                
                error( 'Incomplete specification of mandatory parameters.' );
                
            end
            
        end
        
        function update_pz ( cm )
            
            if ( isempty( cm.b_pz ) )
                
                cm.b_pz = false( cm.alloc_n, 1 );
                cm.pz = zeros( 1, 1, cm.alloc_n );
                
            end
            
            b_todo = cm.b_n & ~cm.b_pz;
            
            n_todo = sum( b_todo );
            
            if ( n_todo > 0 )
                
                cm.update_p_n;
                
                cm.pz( 1, 1, b_todo ) = reshape( sum( cm.p_n( :, b_todo ) .* cm.p_n( :, b_todo ) ), [ 1, 1, n_todo ] );
                
                cm.b_pz = cm.b_n;
                
                if ( cm.len_dp > 0 )
                    
                    tmp = 2 .* reshape( cm.n( b_todo, cm.dp ), [ 1, 1, n_todo, 1, cm.len_dp ] ) .* ...
                        reshape( cm.p_n( :, b_todo )', [ 1, 1, n_todo, 3 ] );
                    
                    % correct for 0 * Inf == NaN
                    
                    tmp( isnan( tmp ) ) = 0;
                    
                    % due to EFn \equiv 0
                    
                    tmp( isinf( tmp ) ) = 0;
                    
                    cm.dpz_dp( 1, 1, b_todo, :, : ) = tmp;
                    
                end
                
            end
            
        end
        
        function update_pxy ( cm, idx_time )
            
            if ( isempty( cm.b_pxy ) )
                
                cm.b_pxy = false( cm.alloc_n, cm.alloc_d );
                cm.pxy = zeros( 2, 1, cm.alloc_n, cm.alloc_d );
                
            end
            
            b_todo = cm.b_n & ~cm.b_pxy( :, idx_time );
            
            n_todo = sum( b_todo );
            
            if ( n_todo > 0 )
                
                cm.update_pz;
                
                tmp = zeros( 2, 1, n_todo );
                
                if ( ~cm.b_s( idx_time ) )
                    
                    if ( cm.b_p( idx_time ) )
                        
                        tmp1 = sum( cm.p( :, idx_time ) .* cm.p( :, idx_time ) ) ./ 3;
                        
                        if ( tmp1 ~= Inf )  %#ok<BDSCI>
                            
                            tmp2 = reshape( sum( cm.p_n( :, b_todo ) .* cm.p( :, idx_time ) ), [ 1, 1, n_todo ] );
                            
                            % correct for Inf * 0 == NaN
                            
                            tmp2( isnan( tmp2 ) ) = 0;
                            
                            tmp( 1, 1, : ) = tmp1 + tmp2;
                            tmp( 2, 1, : ) = tmp1 - tmp2;
                            
                            if ( cm.len_dp > 0 )
                                
                                cm.dpxy_dp( 1, 1, b_todo, :, :, idx_time ) = ...
                                    reshape( cm.n( b_todo, cm.dp ), [ 1, 1, n_todo, 1, cm.len_dp ] ) .* ...
                                    reshape( cm.p( :, idx_time ), [ 1, 1, 1, 3 ] );
                                
                                [ ~, i ] = find( cm.dp == cm.mu( idx_time ) );
                                
                                if ( ~isempty( i ) )
                                    
                                    cm.dpxy_dp( 1, 1, b_todo, :, i, idx_time ) = cm.dpxy_dp( 1, 1, b_todo, :, i, idx_time ) + ...
                                        reshape( cm.p_n( :, b_todo )', [ 1, 1, n_todo, 3 ] );
                                    
                                    cm.dpxy_dp( 2, 1, b_todo, :, :, idx_time ) = - cm.dpxy_dp( 1, 1, b_todo, :, :, idx_time );
                                    
                                    cm.dpxy_dp( :, 1, b_todo, :, i, idx_time ) = cm.dpxy_dp( :, 1, b_todo, :, i, idx_time ) + ...
                                        reshape( 2 .* cm.p( :, idx_time ) ./ 3, [ 1, 1, 1, 3 ] );
                                    
                                else
                                    
                                    cm.dpxy_dp( 2, 1, b_todo, :, :, idx_time ) = - cm.dpxy_dp( 1, 1, b_todo, :, :, idx_time );
                                    
                                end
                                
                            end
                            
                        else
                            
                            tmp( : ) = Inf;
                            
                            if ( cm.len_dp > 0 )
                                
                                cm.dpxy_dp( :, 1, b_todo, :, :, idx_time ) = 0;
                                
                            end
                            
                        end
                        
                    end
                    
                else
                    
                    if ( cm.s( 4, idx_time ) ~= Inf )
                        
                        tmp2 = reshape( sum( cm.p_n( :, b_todo ) .* cm.s( 1 : 3, idx_time ) ), [ 1, 1, n_todo ] );
                        
                        tmp( 1, 1, : ) = cm.s( 4, idx_time ) + tmp2;
                        tmp( 2, 1, : ) = cm.s( 4, idx_time ) - tmp2;
                        
                        if ( cm.len_dp > 0 )
                            
                            cm.dpxy_dp( 1, 1, b_todo, :, :, idx_time ) = ...
                                reshape( cm.n( b_todo, cm.dp ), [ 1, 1, n_todo, 1, cm.len_dp ] ) .* ...
                                reshape( cm.s( 1 : 3, idx_time ), [ 1, 1, 1, 3 ] );
                            
                            cm.dpxy_dp( 2, 1, b_todo, :, :, idx_time ) = - cm.dpxy_dp( 1, 1, b_todo, :, :, idx_time );
                            
                        end
                        
                        if ( cm.len_ds > 0 )
                            
                            [ ~, i ] = find( cm.ds == cm.mu( idx_time ) );
                            
                            if ( ~isempty( i ) )
                                
                                cm.dpxy_ds( 1, 1, b_todo, 1 : 3, i, idx_time ) = ...
                                    reshape( cm.p_n( :, b_todo )', [ 1, 1, n_todo, 3 ] );
                                
                                cm.dpxy_ds( 2, 1, b_todo, 1 : 3, i, idx_time ) = - cm.dpxy_ds( 1, 1, b_todo, 1 : 3, i, idx_time );
                                
                                cm.dpxy_ds( :, 1, b_todo, 4, i, idx_time ) = 1;
                                
                            end
                            
                        end
                        
                    else
                        
                        tmp( : ) = Inf;
                        
                        if ( cm.len_dp > 0 )
                            
                            cm.dpxy_dp( :, 1, b_todo, :, :, idx_time ) = 0;
                            
                        end
                        
                        if ( cm.len_ds > 0 )
                            
                            cm.dpxy_ds( :, 1, b_todo, :, :, idx_time ) = 0;
                            
                        end
                        
                    end
                    
                end
                
                tmp = tmp + cm.pz( 1, 1, b_todo );
                
                % correct for Inf - Inf == NaN
                
                tmp( isnan( tmp ) ) = Inf;
                
                cm.pxy( :, 1, b_todo, idx_time ) = tmp;
                
                % derivatives
                
                if ( cm.len_dp > 0 )
                    
                    cm.dpxy_dp( :, 1, b_todo, :, :, idx_time ) = cm.dpxy_dp( :, 1, b_todo, :, :, idx_time ) + ...
                        cm.dpz_dp( 1, 1, b_todo, :, : );
                    
                end
                
                % we are done
                
                cm.b_pxy( :, idx_time ) = cm.b_n;
                
            end
            
        end
        
        function update_relaxation ( cm, idx_time )
            
            if ( ~cm.b_E( idx_time ) )
                
                cm.E( 1, :, idx_time ) = exp( - cm.R1 .* cm.tau( idx_time ) );
                cm.E( 2, :, idx_time ) = exp( - cm.R2 .* cm.tau( idx_time ) );
                
                % update info
                
                cm.b_E( idx_time ) = true;
                
            end
            
            if ( isempty( cm.b_EFn ) && ...
                    ( max( cm.D ) > 0 || ...
                    cm.len_dD > 0 || ...
                    ( sum( cm.b_p ) > 0 && max( abs( cm.p( : ) ) ) == Inf ) || ...
                    ( sum( cm.b_s ) > 0 && max( abs( cm.s( : ) ) ) == Inf ) ) )
                
                cm.b_EFn = false( cm.alloc_n, cm.alloc_d );
                cm.EFn = zeros( 3, cm.n_tissues, cm.alloc_n, cm.alloc_d );
                
            end
            
            if ( ~isempty( cm.b_EFn ) )
                
                b_todo = cm.b_n & ~cm.b_EFn( :, idx_time );
                
                n_todo = sum( b_todo );
                
                if ( n_todo > 0 )
                    
                    cm.update_pxy( idx_time );
                    
                    tmp = zeros( 3, cm.n_tissues, n_todo );
                    
                    tmp( 1, :, : ) = ...
                        exp( - cm.D .* cm.tau( idx_time ) .* cm.pxy( 1, 1, b_todo, idx_time ) ) .* ...
                        cm.E( 2, :, idx_time );
                    
                    tmp( 2, :, : ) = ...
                        exp( - cm.D .* cm.tau( idx_time ) .* cm.pz( 1, 1, b_todo ) ) .* ...
                        cm.E( 1, :, idx_time );
                    
                    tmp( 3, :, : ) = ...
                        exp( - cm.D .* cm.tau( idx_time ) .* cm.pxy( 2, 1, b_todo, idx_time ) ) .* ...
                        cm.E( 2, :, idx_time );
                    
                    % correct for 0 * Inf == NaN
                    % infinite gradient moments or shapes always refer to spoilers, even for D == 0
                    
                    tmp( isnan( tmp ) ) = 0;
                    
                    cm.EFn( :, :, b_todo, idx_time ) = tmp;
                    
                    % now for the derivatives, if necessary
                    
                    if ( cm.len_dD > 0 )
                        
                        tmp = zeros( 3, cm.dD, n_todo );
                        
                        tmp( 1, :, : ) = ...
                            - cm.tau( idx_time ) .* cm.pxy( 1, 1, b_todo, idx_time ) .* cm.EFn( 1, cm.dD, b_todo, idx_time );
                        
                        tmp( 2, :, : ) = ...
                            - cm.tau( idx_time ) .* cm.pz( 1, 1, b_todo ) .* cm.EFn( 2, cm.dD, b_todo, idx_time );
                        
                        tmp( 3, :, : ) = ...
                            - cm.tau( idx_time ) .* cm.pxy( 2, 1, b_todo, idx_time ) .* cm.EFn( 3, cm.dD, b_todo, idx_time );
                        
                        % resolve Inf * 0 == NaN
                        % The derivative must be zero as EFn \equiv 0 in these cases.
                        
                        tmp( isnan( tmp ) ) = 0;
                        
                        cm.dEFn_dD( :, :, b_todo, idx_time ) = tmp;
                        
                    end
                    
                    if ( cm.len_dtau > 0 )
                        
                        [ ~, i ] = find( cm.dtau == cm.mu( idx_time ) );
                        
                        if ( ~isempty( i ) )
                            
                            tmp = zeros( 3, cm.n_tissues, n_todo );
                            
                            tmp( 1, :, : ) = ...
                                - ( cm.R2 + cm.D .* cm.pxy( 1, 1, b_todo, idx_time ) ) .* cm.EFn( 1, :, b_todo, idx_time );
                            
                            tmp( 2, :, : ) = ...
                                - ( cm.R1 + cm.D .* cm.pz( 1, 1, b_todo ) ) .* cm.EFn( 2, :, b_todo, idx_time );
                            
                            tmp( 3, :, : ) = ...
                                - ( cm.R2 + cm.D .* cm.pxy( 2, 1, b_todo, idx_time ) ) .* cm.EFn( 3, :, b_todo, idx_time );
                            
                            % resolve Inf * 0 == NaN
                            % The derivative must be zero as EFn \equiv 0 in these cases.
                            
                            tmp( isnan( tmp ) ) = 0;
                            
                            cm.dEFn_dtau( :, :, b_todo, i ) = tmp;
                            
                        end
                        
                    end
                    
                    if ( cm.len_dp > 0 )
                        
                        cm.dEFn_dp( 1, :, b_todo, :, :, idx_time ) = ...
                            - cm.tau( idx_time ) .* cm.D .* cm.dpxy_dp( 1, 1, b_todo, :, :, idx_time ) .* cm.EFn( 1, :, b_todo, idx_time );
                        
                        cm.dEFn_dp( 2, :, b_todo, :, :, idx_time ) = ...
                            - cm.tau( idx_time ) .* cm.D .* cm.dpz_dp( 1, 1, b_todo, :, : ) .* cm.EFn( 2, :, b_todo, idx_time );
                        
                        cm.dEFn_dp( 3, :, b_todo, :, :, idx_time ) = ...
                            - cm.tau( idx_time ) .* cm.D .* cm.dpxy_dp( 2, 1, b_todo, :, :, idx_time ) .* cm.EFn( 3, :, b_todo, idx_time );
                        
                    end
                    
                    if ( cm.len_ds > 0 )
                        
                        cm.dEFn_ds( 1, :, b_todo, :, :, idx_time ) = ...
                            - cm.tau( idx_time ) .* cm.D .* cm.dpxy_ds( 1, 1, b_todo, :, :, idx_time ) .* cm.EFn( 1, :, b_todo, idx_time );
                        
                        cm.dEFn_ds( 2, :, b_todo, :, :, idx_time ) = ...
                            - cm.tau( idx_time ) .* cm.D .* cm.dpxy_ds( 2, 1, b_todo, :, :, idx_time ) .* cm.EFn( 3, :, b_todo, idx_time );
                        
                    end
                    
                    % we are done
                    
                    cm.b_EFn( b_todo, idx_time ) = true;
                    
                end
                
                
            end
            
        end
        
        function [ idx_time, b_n_new, b_up_free_time, b_do_free_time ] = pre_time ( cm, handle_time )
            
            if ( cm.options.verbose )
                
                fprintf( 1, 'pre_time: d = %d\n', cm.d );
                fprintf( 1, 'pre_time: n_conf = %d\n', cm.n_conf );
                
            end
            
            % is this a new time interval?
            
            b_time = cm.b_mu & ( cm.mu == handle_time );
            idx_time = find( b_time, 1, 'first' );
            
            if ( isempty( idx_time ) ) % new time index
                
                if ( cm.alloc_d == cm.d ) % need to allocate more space first
                    
                    cm.alloc_d = cm.alloc_d + cm.options.alloc_d;
                    
                    [ ~, rng_d ] = cm.reallocate_nd;
                    
                    % extend b_time
                    
                    b_time( rng_d ) = false;
                    
                end
                
                % get the first free place
                
                idx_time = find( ~cm.b_mu, 1, 'first' );
                
                % update logical indices amd set identifier
                
                b_time( idx_time ) = true;
                cm.b_mu( idx_time ) = true;
                cm.mu( idx_time ) = handle_time;
                
                % configuration order increases by one
                
                cm.d = cm.d + 1;
                
                % prior to the time interval, only the zero order configuration is present
                % and both direct neighbours are targets to be occupied
                
                cm.b_up_free( cm.b_n, idx_time ) = true;
                cm.b_do_free( cm.b_n, idx_time ) = true;
                
            end
            
            % associated with current time interval
            
            b_up_free_time = cm.b_n & cm.b_up_free( :, idx_time );
            b_do_free_time = cm.b_n & cm.b_do_free( :, idx_time );
            
            % how many new states will be created in each direction?
            
            num_up_new = sum( b_up_free_time );
            num_do_new = sum( b_do_free_time );
            num_new = num_up_new + num_do_new;
            
            % determine free space
            
            num_free = sum( ~cm.b_n );
            
            % do we need to allocate more space?
            
            if ( num_free < num_new )
                
                cm.alloc_n = cm.alloc_n + max( num_new - num_free, cm.options.alloc_n );
                
                [ rng_n, ~ ] = cm.reallocate_nd;
                
                % extend b_up_free_time and b_do_free_time
                
                b_up_free_time( rng_n ) = false;
                b_do_free_time( rng_n ) = false;
                
            end
            
            % indices of states which are transferred into previously
            % unoccupied configuration orders
            
            idx_up_old = find( b_up_free_time );
            idx_do_old = find( b_do_free_time );
            
            % corresponding (with respect to order) new indices
            % note that available gaps are filled first
            
            idx_free = find( ~cm.b_n, num_new );
            
            cm.idx_up_new( b_up_free_time ) = idx_free( 1 : num_up_new );
            cm.idx_do_new( b_do_free_time ) = idx_free( num_up_new + 1 : end );
            
            %% bookkeeping for mu == handle_time
            
            % the extremal configurations have new successors
            
            cm.idx_up( b_up_free_time, idx_time ) = cm.idx_up_new( b_up_free_time );
            cm.idx_do( b_do_free_time, idx_time ) = cm.idx_do_new( b_do_free_time );
            
            % same applies (in opposite direction) to the newly created states
            % note that these indices will not be used yet in the current time interval
            
            cm.idx_up( cm.idx_do_new( b_do_free_time ), idx_time ) = idx_do_old;
            cm.idx_do( cm.idx_up_new( b_up_free_time ), idx_time ) = idx_up_old;
            
            % configuration order of handle_time changes
            
            cm.n( cm.idx_up_new( b_up_free_time ), idx_time ) = cm.n( b_up_free_time, idx_time ) + 1;
            cm.n( cm.idx_do_new( b_do_free_time ), idx_time ) = cm.n( b_do_free_time, idx_time ) - 1;
            
            %% bookkeeping for mu ~= handle_time
            
            % associated logical and actual indices
            
            sz_nd = [ cm.alloc_n, cm.alloc_d ];
            
            b_other_mu = cm.b_mu & ~b_time;
            
            b_up_up_o = b_up_free_time & ~cm.b_up_free & b_other_mu;
            b_up_do_o = b_up_free_time & ~cm.b_do_free & b_other_mu;
            b_do_up_o = b_do_free_time & ~cm.b_up_free & b_other_mu;
            b_do_do_o = b_do_free_time & ~cm.b_do_free & b_other_mu;
            
            b_up_up_f = b_up_free_time & cm.b_up_free & b_other_mu;
            b_up_do_f = b_up_free_time & cm.b_do_free & b_other_mu;
            b_do_up_f = b_do_free_time & cm.b_up_free & b_other_mu;
            b_do_do_f = b_do_free_time & cm.b_do_free & b_other_mu;
            
            [ r, c ] = find( b_up_up_o );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                   fprintf( 1, 'b_up_up_o %d x idx_up\n', length( r ) );
                    
                end
                
                ind_up = sub2ind( sz_nd, cm.idx_up_new( r ), c ) ;
                idx_up_other = cm.idx_up( cm.idx_up( b_up_up_o ), idx_time );
                ind_up_other = sub2ind( sz_nd, idx_up_other, c ) ;
                
                cm.idx_up( ind_up ) = idx_up_other;
                cm.idx_do( ind_up_other ) = cm.idx_up_new( r );
                
                cm.b_up_free( ind_up ) = false;
                cm.b_do_free( ind_up_other ) = false;
                
            end
            
            [ r, c ] = find( b_up_do_o );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_up_do_o %d x idx_do\n', length( r ) );
                    
                end

                ind_up = sub2ind( sz_nd, cm.idx_up_new( r ), c ) ;
                idx_up_other = cm.idx_up( cm.idx_do( b_up_do_o ), idx_time );
                ind_up_other = sub2ind( sz_nd, idx_up_other, c ) ;
                
                cm.idx_do( ind_up ) = idx_up_other;
                cm.idx_up( ind_up_other ) = cm.idx_up_new( r );
                
                cm.b_do_free( ind_up ) = false;
                cm.b_up_free( ind_up_other ) = false;                
                
            end
            
            [ r, c ] = find( b_do_up_o );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_do_up_o %d x idx_up\n', length( r ) );
                    
                end
                
                ind_do = sub2ind( sz_nd, cm.idx_do_new( r ), c ) ;
                idx_do_other = cm.idx_do( cm.idx_up( b_do_up_o ), idx_time );
                ind_do_other = sub2ind( sz_nd, idx_do_other, c ) ;
                
                cm.idx_up( ind_do ) = idx_do_other;
                cm.idx_do( ind_do_other ) = cm.idx_do_new( r );
                
                cm.b_up_free( ind_do ) = false;
                cm.b_do_free( ind_do_other ) = false;                

            end
            
            [ r, c ] = find( b_do_do_o );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_do_do_o %d x idx_do\n', length( r ) );
                    
                end
                
                ind_do = sub2ind( sz_nd, cm.idx_do_new( r ), c ) ;
                idx_do_other = cm.idx_do( cm.idx_do( b_do_do_o ), idx_time );
                ind_do_other = sub2ind( sz_nd, idx_do_other, c ) ;
                
                cm.idx_do( ind_do ) = idx_do_other;
                cm.idx_up( ind_do_other ) = cm.idx_do_new( r );
                
                cm.b_do_free( ind_do ) = false;
                cm.b_up_free( ind_do_other ) = false;
                
            end
            
            [ r, c ] = find( b_up_up_f );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_up_up_f %d\n', length( r ) );
                    
                end
                
                cm.b_up_free( sub2ind( sz_nd, cm.idx_up_new( r ), c ) ) = true;
                
            end
            
            [ r, c ] = find( b_up_do_f );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_up_do_f %d\n', length( r ) );
                    
                end
                
                cm.b_do_free( sub2ind( sz_nd, cm.idx_up_new( r ), c ) ) = true;
                
            end
            
            [ r, c ] = find( b_do_up_f );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_do_up_f %d\n', length( r ) );
                    
                end
                
                cm.b_up_free( sub2ind( sz_nd, cm.idx_do_new( r ), c ) ) = true;
                
            end
            
            [ r, c ] = find( b_do_do_f );
            
            if ( ~isempty( r ) )
                
                if ( cm.options.verbose )
                    
%                    fprintf( 1, 'b_do_do_f %d\n', length( r ) );
                    
                end
                
                cm.b_do_free( sub2ind( sz_nd, cm.idx_do_new( r ), c ) ) = true;
                
            end
            
            % configuration order remains unchanged
            
            cm.n( cm.idx_up_new( b_up_free_time ), b_other_mu ) = cm.n( b_up_free_time, b_other_mu );
            cm.n( cm.idx_do_new( b_do_free_time ), b_other_mu ) = cm.n( b_do_free_time, b_other_mu );
            
            %% logical indices of new states
            
            cm.b_up_free( cm.idx_up_new( b_up_free_time ), idx_time ) = true;
            cm.b_do_free( cm.idx_up_new( b_up_free_time ), idx_time ) = false;
            
            cm.b_do_free( cm.idx_do_new( b_do_free_time ), idx_time ) = true;
            cm.b_up_free( cm.idx_do_new( b_do_free_time ), idx_time ) = false;
            
            %% the following is a bit subtle...
            
%             b_up_occ_time = cm.b_n & ~b_up_free_time;
%             b_do_occ_time = cm.b_n & ~b_do_free_time;
%             
%             b_up_free_other = b_other_mu & cm.b_up_free;
%             b_do_free_other = b_other_mu & cm.b_do_free;
%             
%             % first, we look at the cases, for which we can easily determine,
%             % whether the free neighbor will be occupied
%             
%             [ r, c ] = find( b_up_occ_time & b_up_free_other );
%             ind = sub2ind( sz_nd, r, c );
%             r_up = cm.idx_up( r, idx_time );
%             ind_up = sub2ind( sz_nd, r_up, c );
%             b_tmp = ~cm.b_up_free( ind_up );
%             
%             if ( cm.options.verbose )
%                 
%                 fprintf( 1, 'up_up: length( r ) = %d, sum( b_tmp ) = %d\n', length( r ), sum( b_tmp ) );
%                 
%             end
%             
%             cm.b_up_free( ind( b_tmp ) ) = false;
%             r_new = cm.idx_do_new( cm.idx_up( ind_up( b_tmp) ) );
%             cm.idx_up( ind( b_tmp) ) = r_new;
%             cm.idx_do( sub2ind( sz_nd, r_new, c( b_tmp ) ) ) = r( b_tmp );
%             
%             [ r, c ] = find( b_up_occ_time & b_do_free_other );
%             ind = sub2ind( sz_nd, r, c );
%             r_up = cm.idx_up( r, idx_time );
%             ind_up = sub2ind( sz_nd, r_up, c );
%             b_tmp = ~cm.b_do_free( ind_up );
% 
%             if ( cm.options.verbose )
%                 
%                 fprintf( 1, 'up_do: length( r ) = %d, sum( b_tmp ) = %d\n', length( r ), sum( b_tmp ) );
%                 
%             end
%             
%             cm.b_do_free( ind( b_tmp ) ) = false;
%             r_new = cm.idx_do_new( cm.idx_do( ind_up( b_tmp) ) );
%             cm.idx_do( ind( b_tmp) ) = r_new;
%             cm.idx_up( sub2ind( sz_nd, r_new, c( b_tmp ) ) ) = r( b_tmp );
%                         
%             [ r, c ] = find( b_do_occ_time & b_up_free_other );
%             ind = sub2ind( sz_nd, r, c );
%             r_do = cm.idx_do( r, idx_time );
%             ind_do = sub2ind( sz_nd, r_do, c );
%             b_tmp = ~cm.b_up_free( ind_do );
%             
%             if ( cm.options.verbose )
%                 
%                 fprintf( 1, 'do_up: length( r ) = %d, sum( b_tmp ) = %d\n', length( r ), sum( b_tmp ) );
%                 
%             end
%             
%             cm.b_up_free( ind( b_tmp ) ) = false;
%             r_new = cm.idx_up_new( cm.idx_up( ind_do( b_tmp) ) );
%             cm.idx_up( ind( b_tmp) ) = r_new;
%             cm.idx_do( sub2ind( sz_nd, r_new, c( b_tmp ) ) ) = r( b_tmp );
%             
%             [ r, c ] = find( b_do_occ_time & b_do_free_other );
%             ind = sub2ind( sz_nd, r, c );
%             r_do = cm.idx_do( r, idx_time );
%             ind_do = sub2ind( sz_nd, r_do, c );
%             b_tmp = ~cm.b_do_free( ind_do );
% 
%             if ( cm.options.verbose )
%                 
%                 fprintf( 1, 'do_do: length( r ) = %d, sum( b_tmp ) = %d\n', length( r ), sum( b_tmp ) );
%                 
%             end
%             
%             cm.b_do_free( ind( b_tmp ) ) = false;
%             r_new = cm.idx_up_new( cm.idx_do( ind_do( b_tmp) ) );
%             cm.idx_do( ind( b_tmp) ) = r_new;
%             cm.idx_up( sub2ind( sz_nd, r_new, c( b_tmp ) ) ) = r( b_tmp );
            
            % now we can update this as well
            
            cm.b_up_free( b_up_free_time, idx_time ) = false;
            cm.b_do_free( b_do_free_time, idx_time ) = false;
            
            %% prepare to update the new occupation
            
            b_n_new = cm.b_n;
            b_n_new( cm.idx_up_new( b_up_free_time ) ) = true;
            b_n_new( cm.idx_do_new( b_do_free_time ) ) = true;
            
        end
        
        function post_time ( cm, b_n_new )
            
            cm.b_n = b_n_new;
            cm.n_conf = sum( cm.b_n );
            
            if ( cm.options.verbose )
                
                fprintf( 1, 'post_time: d = %d\n', cm.d );
                fprintf( 1, 'post_time: n_conf = %d\n', cm.n_conf );
                
            end
            
        end
                
        function b_is_ok = meltdown ( cm )
            %% remove negligible configurations
            
            b_is_ok = true;
            
            abs_m_2 = zeros( cm.alloc_n, 1 );
            abs_m_2( cm.b_n ) = real( sum( reshape( real( cm.m( :, :, cm.b_n ) .* conj( cm.m( :, :, cm.b_n ) ) ), [ 3 * cm.n_tissues, cm.n_conf ] ) ) ).';
            b_small = abs_m_2 < cm.epsilon^2;

            % needed a lot below
            
            sz_nd = [ cm.alloc_n, cm.alloc_d ];
            
            % we do not want to remove these (to be updated later)
            
            b_exclude_n = ~cm.b_n;
            b_exclude_n( cm.null_idx ) = true;
            
            while( true )
                
                %% candidates for removal:
                % - occupied
                % - not explicitly excluded
                % - small
                % - have an unoccupied neighbor in every direction
                
                b_remove_n = ...
                    cm.b_n & ...
                    ~b_exclude_n & ...
                    b_small & ...
                    ( sum( cm.b_up_free( :, cm.b_mu ) | cm.b_do_free( :, cm.b_mu ), 2 ) == cm.d );
                
                n_remove = sum( b_remove_n );
                
                if ( n_remove == 0 )
                    
                    break;  % nothing more to to
                    
                end
                
                removed = 0;
                
                if ( cm.rapid_meltdown ) % looking at the candidates at once (requires more memory)
                
                    % first we have to check that no occupation patterns like
                    %          ---------
                    %          | 1 | 0 |
                    %          ---------
                    %          | 0 | 1 |
                    %          ---------
                    % arise after removal of the configurations

                    idx_n = find( b_remove_n );
                    idx_mu = find( cm.b_mu );

                    b_up_occ = ~cm.b_up_free( b_remove_n, cm.b_mu );
                    b_do_occ = ~cm.b_do_free( b_remove_n, cm.b_mu );

                    idx_up_occ = cm.idx_up( b_remove_n, cm.b_mu );
                    idx_do_occ = cm.idx_do( b_remove_n, cm.b_mu );

                    % up up

                    mu_pairs = ...
                        reshape( idx_mu, [ 1, cm.d ] ) > reshape( idx_mu, [ 1, 1, cm.d ] );
                    
                    [ r_up_up, c_up_up ] = ...
                        find( reshape( ...
                            b_up_occ & reshape( b_up_occ, [ n_remove, 1, cm.d ] ) & mu_pairs, ...
                            [ n_remove, cm.d * cm.d ] ) );

                    [ i_up, j_up ] = ind2sub( [ cm.d, cm.d ], c_up_up );

                    n_ = length( i_up );
                    sz_n_ = [ n_, 1 ];
                    
                    if ( n_ > 0 )
                        
                        b_free = cm.b_up_free( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_up_occ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_up_up( : ), ...
                            i_up( : ) ) ), sz_n_ ), ...
                            reshape( idx_mu( j_up ), sz_n_ ) ) );
                        
                        b_exclude_n( idx_n( r_up_up( b_free ) ) ) = true;
                        
                        sz_n_ = [ sum( ~b_free ), 1 ];
                        
                        b_exclude_n( cm.idx_up( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_up_occ ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_up_up( ~b_free), ...
                            i_up( ~b_free ) ) ), sz_n_ ), ...
                            reshape( idx_mu( j_up( ~b_free ) ), sz_n_ ) ) ) ) = true;
                        
                    end
                    
                    % do do
                    
                    [ r_do_do, c_do_do ] = ...
                        find( reshape( ...
                            b_do_occ & reshape( b_do_occ, [ n_remove, 1, cm.d ] ) & mu_pairs, ...
                            [ n_remove, cm.d * cm.d ] ) );

                    [ i_do, j_do ] = ind2sub( [ cm.d, cm.d ], c_do_do );

                    n_ = length( i_do );
                    sz_n_ = [ n_, 1 ];

                    if ( n_ > 0 )
                        
                        b_free = cm.b_do_free( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_do_occ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_do_do( : ), ...
                            i_do( : ) ) ), sz_n_ ) , ...
                            reshape( idx_mu( j_do ), sz_n_ )  ) );
                        
                        b_exclude_n( idx_n( r_do_do( b_free ) ) ) = true;

                        sz_n_ = [ sum( ~b_free ), 1 ];
                                                
                        b_exclude_n( cm.idx_do( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_do_occ ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_do_do( ~b_free), ...
                            i_do( ~b_free ) ) ), sz_n_ ) , ...
                            reshape( idx_mu( j_do( ~b_free ) ), sz_n_ ) ) ) ) = true;
                        
                    end
                             
                    % up do
                    
                    mu_pairs = ...
                        reshape( idx_mu, [ 1, cm.d ] ) ~= reshape( idx_mu, [ 1, 1, cm.d ] );
                    
                    [ r_up_do, c_up_do ] = ...
                        find( reshape( ...
                            b_up_occ & reshape( b_do_occ, [ n_remove, 1, cm.d ] ) & mu_pairs, ...
                            [ n_remove, cm.d * cm.d ] ) );

                    [ i_up, j_do ] = ind2sub( [ cm.d, cm.d ], c_up_do );

                    n_ = length( i_up );
                    sz_n_ = [ n_, 1 ];
                                        
                    if ( n_ > 0 )
                        
                        b_free = cm.b_do_free( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_up_occ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_up_do( : ), ...
                            i_up( : ) ) ), sz_n_ ), ...
                            reshape( idx_mu( j_do ), sz_n_ )  ) );
                        
                        b_exclude_n( idx_n( r_up_do( b_free ) ) ) = true;
                        
                        sz_n_ = [ sum( ~b_free ), 1 ];
                        
                        b_exclude_n( cm.idx_do( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_up_occ ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_up_do( ~b_free), ...
                            i_up( ~b_free ) ) ), sz_n_ ) , ...
                            reshape( idx_mu( j_do( ~b_free ) ), sz_n_ ) ) ) ) = true;
                        
                    end
                    
                    % update configurations to be removed
                    
                    b_remove_n( b_exclude_n ) = false;
                    removed = sum( b_remove_n );
                    
                    % now we may remove the configurations
                    
                    cm.b_n( b_remove_n ) = false;
                    cm.n_conf = cm.n_conf - removed;
                    
                    % disconnect adjacent configurations
                    
                    b_up_occ = b_remove_n & cm.b_mu & ~cm.b_up_free;
                    b_do_occ = b_remove_n & cm.b_mu & ~cm.b_do_free;
                    
                    [ r_up_occ, c_up_occ ] = find( b_up_occ );
                    [ r_do_occ, c_do_occ ] = find( b_do_occ );
                    
                    idx_up_occ = cm.idx_up( b_up_occ );
                    idx_do_occ = cm.idx_do( b_do_occ );
                    
                    cm.b_up_free( sub2ind( sz_nd, idx_do_occ, c_do_occ ) ) = true;
                    cm.b_do_free( sub2ind( sz_nd, idx_up_occ, c_up_occ ) ) = true;
                    cm.b_up_free( sub2ind( sz_nd, r_up_occ, c_up_occ ) ) = true;
                    cm.b_do_free( sub2ind( sz_nd, r_do_occ, c_do_occ ) ) = true;
                    
                    if ( ~isempty( cm.b_tau_n ) )
                        
                        cm.b_tau_n( b_remove_n ) = false;
                        
                    end
                    
                    if ( ~isempty( cm.b_p_n ) )
                        
                        cm.b_p_n( b_remove_n ) = false;
                        
                    end
                    
                    if ( ~isempty( cm.b_pz ) )
                        
                        cm.b_pz( b_remove_n ) = false;
                        
                    end
                    
                    if ( ~isempty( cm.b_pxy ) )
                        
                        cm.b_pxy( b_remove_n, cm.b_mu ) = false;
                        
                    end
                    
                    if ( ~isempty( cm.b_EFn ) )
                        
                        cm.b_EFn( b_remove_n, cm.b_mu ) = false;
                        
                    end

                else % look at the candidates individually
                    
                    r = find( b_remove_n );
                
                    for i = 1 : n_remove
                        
                        % first we have to check that no occupation patterns like
                        %          ---------
                        %          | 1 | 0 |
                        %          ---------
                        %          | 0 | 1 |
                        %          ---------
                        % arise after removal of the actual configuration r( i )
                        
                        b_up_occ = cm.b_mu & ~cm.b_up_free( r( i ), : );
                        b_do_occ = cm.b_mu & ~cm.b_do_free( r( i ), : );
                        
                        c_up_occ = find( b_up_occ );
                        c_do_occ = find( b_do_occ );
                        
                        n_up_occ = length( c_up_occ );
                        n_do_occ = length( c_do_occ );
                        
                        idx_up_occ = cm.idx_up( r( i ), c_up_occ )';
                        idx_do_occ = cm.idx_do( r( i ), c_do_occ )';
                        
                        c_up_occ = c_up_occ( : );
                        c_do_occ = c_do_occ( : );
                        
                        [ r_, c_ ] = find( ( 1 : n_up_occ ) > ( 1 : n_up_occ)' );
                        
                        if ( any( cm.b_up_free( sub2ind( sz_nd, idx_up_occ( r_( : ) ), c_up_occ( c_( : ) ) ) ) ) )
                            
                            b_exclude_n( r( i ) ) = true;
                            continue;
                            
                        end
                        
                        [ r_, c_ ] = find( true( n_up_occ, n_do_occ ) );
                        
                        if ( any( cm.b_up_free( sub2ind( sz_nd, idx_do_occ( c_( : ) ), c_up_occ( r_( : ) ) ) ) ) )
                            
                            b_exclude_n( r( i ) ) = true;
                            continue;
                            
                        end
                        
                        [ r_, c_ ] = find( ( 1 : n_do_occ ) > ( 1 : n_do_occ)' );
                        
                        if ( any( cm.b_do_free( sub2ind( sz_nd, idx_do_occ( r_ ), c_do_occ( c_ ) ) ) ) )
                            
                            b_exclude_n( r( i ) ) = true;
                            continue;
                            
                        end
                        
                        % now we may remove the configuration r( i )
                        
                        removed = removed + 1;
                        
                        cm.b_n( r( i ) ) = false;
                        cm.n_conf = cm.n_conf - 1;
                        
                        % disconnect with adjacent configurations
                        
                        cm.b_up_free( sub2ind( sz_nd, idx_do_occ, c_do_occ ) ) = true;
                        cm.b_do_free( sub2ind( sz_nd, idx_up_occ, c_up_occ ) ) = true;
                        cm.b_up_free( r( i ), cm.b_mu ) = true;
                        cm.b_do_free( r( i ), cm.b_mu ) = true;
                        
                        if ( ~isempty( cm.b_tau_n ) )
                            
                            cm.b_tau_n( r( i ) ) = false;
                            
                        end
                        
                        if ( ~isempty( cm.b_p_n ) )
                            
                            cm.b_p_n( r( i ) ) = false;
                            
                        end
                        
                        if ( ~isempty( cm.b_pz ) )
                            
                            cm.b_pz( r( i ) ) = false;
                            
                        end
                        
                        if ( ~isempty( cm.b_pxy ) )
                            
                            cm.b_pxy( r( i ), cm.b_mu ) = false;
                            
                        end
                        
                        if ( ~isempty( cm.b_EFn ) )
                            
                            cm.b_EFn( r( i ), cm.b_mu ) = false;
                            
                        end
                                             
                    end
                    
                end
                
                % dimensions without nonzero configurations are in fact
                % nonexistent and can be removed
                
                b_remove_mu = cm.b_mu & ~any( cm.n( cm.b_n, : ) ~= 0 );
                
                if ( any( b_remove_mu ) )
                    
                    if ( cm.options.verbose )
                        
                        fprintf( 1, 'meltdown: Remove %d dimensions.\n', sum( b_remove_mu ) );
                        
                    end
                    
                    cm.b_mu( b_remove_mu ) = false;
                    cm.d = sum( cm.b_mu );
                    
                    cm.b_up_free( :, b_remove_mu ) = true;
                    cm.b_do_free( :, b_remove_mu ) = true;
                    
                    cm.b_tau( b_remove_mu ) = false;
                    cm.b_p( b_remove_mu ) = false;
                    cm.b_s( b_remove_mu ) = false;
                    
                    cm.b_E( b_remove_mu ) = false;
                    
                    if ( ~isempty( cm.b_pxy ) )
                        
                        cm.b_pxy( :, b_remove_mu ) = false;
                        
                    end
                    
                    if ( ~isempty( cm.b_EFn ) )
                        
                        cm.b_EFn( :, b_remove_mu ) = false;
                        
                    end
                    
                end                
                
                if ( cm.options.debug )
                    
                    str = [ 'after removing ', num2str( removed ), ' configurations' ];
                    
                    if ( ~cm.check_state( str ) )
                        
                        b_is_ok = false;
                        return;
                        
                    end
                    
                end
                
            end
            
        end
        
        function [ rng_n, rng_d ] = reallocate_nd ( cm )
            
            [ alloc_n_old, alloc_d_old ] = size( cm.n );
            
            if ( cm.options.verbose )
                
                fprintf( 1, 'reallocate_nd: (before) alloc_n = %d\n', alloc_n_old );
                fprintf( 1, 'reallocate_nd: (before) alloc_d = %d\n', alloc_d_old );
                
            end
            
            rng_n = alloc_n_old + 1 : cm.alloc_n;
            rng_d = alloc_d_old + 1 : cm.alloc_d;
            
            if ( isempty( rng_n ) && isempty( rng_d ) )
                
                error( 'No reason to reallocate.' );
                
            end
            
            if ( ~isempty( rng_n ) )
                
                cm.b_n( rng_n ) = false;
                cm.n( rng_n, : ) = 0;
                
                cm.m( :, :, rng_n ) = 0;
                
                cm.b_up_free( rng_n, : ) = true;
                cm.b_do_free( rng_n, : ) = true;
                
                cm.idx_up( rng_n, : ) = 0;
                cm.idx_do( rng_n, : ) = 0;
                
                cm.idx_up_new( rng_n ) = 0;
                cm.idx_do_new( rng_n ) = 0;
                
                if ( ~isempty( cm.b_tau_n ) )
                    
                    cm.b_tau_n( rng_n ) = false;
                    cm.tau_n( rng_n ) = 0;
                    
                end
                
                if ( ~isempty( cm.b_p_n ) )
                    
                    cm.b_p_n( rng_n ) = false;
                    cm.p_n( :, rng_n ) = 0;
                    
                end
                
                if ( ~isempty( cm.b_pz ) )
                    
                    cm.b_pz( rng_n ) = false;
                    cm.pz( :, :, rng_n ) = 0;
                    
                end
                
                if ( ~isempty( cm.b_pxy ) )
                    
                    cm.b_pxy( rng_n, : ) = false;
                    cm.pxy( :, :, rng_n, : ) = 0;
                    
                end
                
                if ( ~isempty( cm.b_EFn ) )
                    
                    cm.b_EFn( rng_n, : ) = false;
                    cm.EFn( rng_n, : ) = 0;
                    
                end
                
                if ( cm.len_dR1 > 0 )
                    
                    cm.dm_dR1( :, :, rng_n ) = 0;
                    
                end
                
                if ( cm.len_dR2 > 0 )
                    
                    cm.dm_dR2( :, :, rng_n ) = 0;
                    
                end
                
                if ( cm.len_dD > 0 )
                    
                    cm.dm_dD( :, :, rng_n ) = 0;
                    cm.dEFn_dD( :, :, rng_n, : ) = 0;
                    
                end
                
                if ( cm.len_dB1 > 0 )
                    
                    cm.dm_dB1( :, :, rng_n ) = 0;
                    
                end
                
                if ( cm.len_dFlipAngle > 0 )
                    
                    cm.dm_dFlipAngle( :, :, rng_n, : ) = 0;
                    
                end
                
                if ( cm.len_dPhase > 0 )
                    
                    cm.dm_dPhase( :, :, rng_n, : ) = 0;
                    
                end
                
                if ( cm.len_dtau > 0 )
                    
                    cm.dm_dtau( :, :, rng_n, : ) = 0;
                    cm.dEFn_dtau( :, :, rng_n, : ) = 0;
                    
                end
                
                if ( cm.len_dp > 0 )
                    
                    cm.dm_dp( :, :, rng_n, :, : ) = 0;
                    cm.dpz_dp( :, :, rng_n, :, : ) = 0;
                    cm.dpxy_dp( :, :, rng_n, :, :, : ) = 0;
                    cm.dEFn_dp( :, :, rng_n, :, :, : ) = 0;
                    
                end
                
                if ( cm.len_ds > 0 )
                    
                    cm.dm_ds( :, :, rng_n, :, : ) = 0;
                    cm.dpxy_ds( :, :, rng_n, :, :, : ) = 0;
                    cm.dEFn_ds( :, :, rng_n, :, :, : ) = 0;
                    
                end
                
            end
            
            if ( ~isempty( rng_d ) )
                
                cm.b_mu( rng_d ) = false;
                cm.mu( rng_d ) = 0;
                
                cm.n( :, rng_d ) = 0;
                
                cm.b_up_free( :, rng_d ) = true;
                cm.b_do_free( :, rng_d ) = true;
                
                cm.idx_up( :, rng_d ) = 0;
                cm.idx_do( :, rng_d ) = 0;
                
                cm.b_tau( rng_d ) = false;
                cm.tau( rng_d, 1 ) = 0;
                
                cm.b_p( rng_d ) = false;
                cm.p( :, rng_d ) = 0;
                
                cm.b_s( rng_d ) = false;
                cm.s( :, rng_d ) = 0;
                
                cm.b_E( rng_d ) = false;
                cm.E( :, :, rng_d ) = 0;
                
                if ( ~isempty( cm.b_pxy ) )
                    
                    cm.b_pxy( :, rng_d ) = false;
                    cm.pxy( :, :, :, rng_d ) = 0;
                    
                end
                
                if ( ~isempty( cm.b_EFn ) )
                    
                    cm.b_EFn( :, rng_d ) = false;
                    cm.EFn( :, rng_d ) = 0;
                    
                end
                
                if ( cm.len_dD > 0 )
                    
                    cm.dEFn_dD( :, :, :, rng_d ) = 0;
                    
                end
                
                if ( cm.len_dp > 0 )
                    
                    cm.dpxy_dp( :, :, :, :, :, rng_d ) = 0;
                    cm.dEFn_dp( :, :, :, :, :, rng_d ) = 0;
                    
                end
                
                if ( cm.len_ds > 0 )
                    
                    cm.dpxy_ds( :, :, :, :, :, rng_d ) = 0;
                    cm.dEFn_ds( :, :, :, :, :, rng_d ) = 0;
                    
                end
                
            end
            
            if ( cm.options.verbose )
                
                fprintf( 1, 'reallocate_nd: (after) alloc_n = %d\n', cm.alloc_n );
                fprintf( 1, 'reallocate_nd: (after) alloc_d = %d\n', cm.alloc_d );
                
            end
            
        end
        
    end
    
end