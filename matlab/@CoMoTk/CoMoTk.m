classdef CoMoTk < matlab.mixin.Copyable
    % CoMoTk  Configuration Model Toolkit
    %
    % Implementation of the configuration model (CM), which has
    % been extensively described in a separate open-access manuscript.
    % (link available on: https://github.com/cganter/CoMoTk)
    % 
    % Also check out the specific documentation and example scripts.
    %
    % To facilitate the use in optimization problems, it also supports the
    % (optional) calculation of several first order derivatives.
    %
    % Requires Matlab R2016b (or later).
    %
    % Carl Ganter 2018
    
    properties ( Constant )
        
        % frequently used parameters
        
        sqrt_2 = sqrt( 2 );
        sqrt_0p5 = sqrt( 0.5 );
        
        % diffusion tensor indices

        Tensor_idx = struct( ...
			'xx',      1, ...
			'yy',      2, ...
			'zz',      3, ...
			'xy',      4, ...
			'xz',      5, ...
			'yz',      6 ...
			);

        % diffusion variants
        
        No_Diff = 1;
        Iso_Diff = 2;
        Tensor_Diff = 3;
        
    end
    
    properties
        
        % function handle for inhomogeneous broadening
        %
        % the format is 'inhomogeneous_decay( t, param )'
        % where 't' is the time and 'param' may contain arbitrary fields with
        % optional parameters (supplied in the equally named argument of
        % the method 'sum')
        % Any user supplied function needs to accept these two parameters
        % (if 'param' is not needed, its value can simply be ignored)
        
        inhomogeneous_decay = [];
        
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
        
        mu;
        
        % delta omega (due to chemical shift, e.g. to specify water/fat models)
        
        dom;                                                               % [ rad / ms ]
        
        % transition rates (magnetization transfer and exchange)
        
        k;                                                                 % [ 1 / ms ]
        
        % as unit for bulk motion we assume [ mm / s ]
        % (used in method "time")
        
        %% measurement conditions
        % relative B1
        
        B1;
                
    end
    
    properties ( SetAccess = private )
        
        b_prep = false;                                                    % is everything prepared?
        
        %% global options
        % the purpose of the alloc_* variables is to increase the numerical performance by minimizing
        % -> recalculation of identical quantities (alloc_RF)
        % -> reallocation of memory (alloc_d, alloc_n)
        % the set values define
        % -> the initial state
        % -> and the minimal reallocation increment, if a limit is reached
        
        options_priv = struct( ...
            'alloc_d', 1, ...                                              % dimension of configuration model
            'alloc_n', 10000, ...                                          % number of occupied states
            'epsilon', 0, ...                                              % discard configuration vectors with a squared norm smaller than this
            'rapid_meltdown', true, ...                                    % should be faster, but requires more memory
            'verbose', false, ...                                          % various status messages
            'debug', false ...                                             % store state for debugging purposes
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
        
        b_lambda = [];                                                     % occupied dimensions
        b_n = [];                                                          % occupied configurations (total)
        
        lambda = [];                                                       % indices of dimensions
        n = [];                                                            % configuration orders in each dimension
        
        n_tissues = 0;		                                               % number of tissues
        d = 0;                                                             % == sum( cm.b_lambda ) : dimension of configuration model
        n_conf = 0;                                                        % == sum( cm.b_n ) : occupied configurations (including n == 0)
        null_idx = 1;                                                      % index, corresponding to n == 0 (remains untouched)
        
        %% mandatory tissue properties
        % general properties
        
        diffusion = CoMoTk.No_Diff;
        
        % cf. dependent counterparts
        
        R1_priv = [];
        R2_priv = [];
        D_priv = [];
           
        %% optional tissue properties and measurement conditions
        % cf. dependent counterparts
        
        mu_priv = 1;
        dom_priv = 0;
        k_priv = [];
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
        p = [];                                                            % [ rad / um ]
        
        % integrals for arbitrary (but otherwise constant) gradient shapes
        % see documentation for the precise definition
        
        s = [];                                                            % [ rad / um ]
        S = [];                                                            % [ rad^2 / um^2 ]
   
        % auxiliary variable for tensor diffusion
        
        D_vec = [];
        
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
        dS = [];
        
        % associated numbers ...
                
        len = struct( ...
            'dR1', 0, ...
            'dR2', 0, ...
            'dD', 0, ...
            'dB1', 0, ...
            'dFlipAngle', 0, ...
            'dPhase', 0, ...
            'dtau', 0, ...
            'dp', 0, ...
            'ds', 0, ...
            'dS', 0 ...
            );
        
        % ... degrees of freedom ...
        
        dim = struct( ...
            'R1', 1, ...
            'R2', 1, ...
            'D', 1, ...
            'B1', 1, ...
            'FlipAngle', 1, ...
            'Phase', 1, ...
            'tau', 1, ...
            'p', 3, ...
            's', 3, ...
            'S', 1 ...
        );
        
        % ... and derivatives of configuration vector
        
        dm = struct( ...
            'dR1', [], ...
            'dR2', [], ...
            'dD', [], ...
            'dB1', [], ...
            'dFlipAngle', [], ...
            'dPhase', [], ...
            'dtau', [], ...
            'dp', [], ...
            'ds', [], ...
            'dS', [] ...
            );
        
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
        
        % relaxation with diffusion and derivatives
        
        b_En = [];
        En = [];
        En_exp = [];
        
        dEn_exp_dD = [];
        dEn_exp_dtau = [];
        dEn_exp_dp = [];
        dEn_exp_ds = [];
        dEn_exp_dS = [];
        
        % debugging
        
        log = [];
        
    end
    
    methods
        
        %% constructor
        
        function cm = CoMoTk ( )
            % constructor (no arguments, no action)
            
        end
        
        %% set / get methods (mandatory variables)
        
        function set.R1 ( cm, R1 )
            % set relaxation rate R1 = 1 / T1
            %
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( min( R1 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R1 );
                cm.mu = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( length( R1 ) ~= cm.n_tissues )
                
                error( 'length( R1 ) must be equal to number of tissues.' );
                
            end
            
            cm.R1_priv = reshape( R1, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.R1 ( cm )
            % return relaxation rate
            
            res = cm.R1_priv;
            
        end
        
        function set.R2 ( cm, R2 )
            % set relaxation rate R2 = 1 / T2
            %
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( min( R2 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R2 );
                cm.mu = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( length( R2 ) ~= cm.n_tissues )
                
                error( 'length( R2 ) must be equal to number of tissues.' );
                
            end
            
            cm.R2_priv = reshape( R2, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.R2 ( cm )
            % return relaxation rate
            
            res = cm.R2_priv;
            
        end
        
        function set.D ( cm, D )
            % set diffusion constant or tensor
            %
            % determined by size of argument:
            % size( D ) == [ 1, n_tissues ]    (isotropic diffusion constant)
            % size( D ) == [ 3, 3, n_tissues ] (anisotropic diffusion tensor)
            %
            % absence of diffusion: just set all elements D( : ) to zero 
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            sz_D = size( D );
            
            if ( sz_D( 1 ) ~= 3 || sz_D( 2 ) ~= 3 )
                
                cm.diffusion = CoMoTk.Iso_Diff;
                
                n_tis = sz_D( 2 );
                
            else
                
                cm.diffusion = CoMoTk.Tensor_Diff;
                
                if ( length( sz_D ) == 3 )
                    
                    n_tis = sz_D( 3 );
                    
                else
                    
                    n_tis = 1;
                    
                end
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = n_tis;
                cm.mu = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( n_tis ~= cm.n_tissues )
                
                error( 'Supplied D does not match number of tissues.' );
                
            end
            
            if ( cm.diffusion == CoMoTk.Tensor_Diff )
                
                cm.D_vec = zeros( 6, cm.n_tissues );
                
                cm.D_vec( CoMoTk.Tensor_idx.xx, : ) = D( 1, 1, : );
                cm.D_vec( CoMoTk.Tensor_idx.yy, : ) = D( 2, 2, : );
                cm.D_vec( CoMoTk.Tensor_idx.zz, : ) = D( 3, 3, : );
                cm.D_vec( CoMoTk.Tensor_idx.xy, : ) = 2 .* D( 1, 2, : );
                cm.D_vec( CoMoTk.Tensor_idx.xz, : ) = 2 .* D( 1, 3, : );
                cm.D_vec( CoMoTk.Tensor_idx.yz, : ) = 2 .* D( 2, 3, : );
                
            end
            
            cm.D_priv = D;
            
            if ( max( abs( cm.D( : ) ) ) == 0 )
                
                cm.diffusion = CoMoTk.No_Diff;
                
            end
            
        end
        
        function res = get.D ( cm )
            % return diffusion constant or tensor
            
            res = cm.D_priv;
            
        end
        
        %% set / get methods (optional variables)
        
        function set.B1 ( cm, B1 )
            % set relative B1+
            %
            % set to 1, if not explicitly specified            
            
            if ( cm.b_prep )
                
                error( 'B1 can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( B1 < 0 )
                
                error( 'B1 must not be negative.' );
                
            end
            
            cm.B1_priv = B1;
            
        end
        
        function res = get.B1 ( cm )
            % get relative B1+
            
            res = cm.B1_priv;
            
        end
        
        % (relative) proton density of complex tissues
        % (default == 1)
        
        function set.mu ( cm, mu )
            % set the (relative) proton density
            %
            % default = 1 (only, if cm.n_tissues == 1)
            % mandatory parameter, if cm.n_tissues > 1
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMoTk.prep.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( mu );
                
            elseif ( length( mu ) ~= cm.n_tissues )
                
                error( 'length( mu ) must be equal to number of tissues.' );
                
            end
            
            cm.mu_priv = reshape( mu, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.mu ( cm )
            % return the (relative) proton density
            
            res = cm.mu_priv;
            
        end
        
        function set.dom ( cm, dom )
            % relative resonance frequencies of subtissues 
            % (e.g. due to chemical shift)
            %
            % length must match the number of subtissues (e.g. for fat models)            
            % set to 0, if not explicitly specified
                        
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( dom );
                
            elseif ( length( dom ) ~= cm.n_tissues )
                
                error( 'length( dom ) must be equal to number of tissues.' );
                
            end
            
            cm.dom_priv = reshape( dom, [ cm.n_tissues, 1 ] );
            
        end
        
        function res = get.dom ( cm )
            % return relative resonance frequencies
            
            res = cm.dom_priv;
            
        end
                
        function set.k ( cm, k )
            % transition rates for magnetization transfer/exchange
            %
            % the following size is required:
            % size( k ) == [ n_tissues, n_tissues ]
            % set to [], if not set (non-exchanging tissues)
            % 
            % Note: 
            % Diagonal elements are calculated according to sum rule
            % (particle conservation). Supplied diagonal values are ignored. 
            
            if ( cm.n_tissues == 0 )
                
                error( 'initialize tissues before setting the transition rates');
                
            end 
            
            sz_k = size( k );
            
            if ( length( sz_k ) ~= 2 || sz_k( 1 ) ~= cm.n_tissues || sz_k( 2 ) ~= cm.n_tissues )
                
                error( 'size( k ) ~= [ n_tissues, n_tissues ].' );
                
            end
                        
            % establish sum rule (particle conservation)

            k = k - diag( diag( k ) );  % first, set the diagonal to zero

            cm.k_priv = k - diag( sum( k ) );   % now we assign the negative row sum to the diagonal
                                                % and store the result
                            
        end
        
        function res = get.k( cm )
            % get transition rates (== [], if not set)
            
            res = cm.k_priv;
            
        end
        
        % options
        
        function set.options ( cm, options )
            % set some options
            
            if ( cm.b_prep )
                
                error( 'Too late to change options.' );
                
            end
            
            cm.options_priv = options;
            
        end
        
        function res = get.options ( cm )
            % get actual options
            %
            % can be used to obtain actual defaults, before applying the
            % associated 'set' method
            
            res = cm.options_priv;
            
        end
        
        %% define, which derivatives should be calculated
        
        function set_derivatives ( cm, param )
            % Define the set of derivatives, to be calculated
            %
            % IN
            %
            % every parameter X, for which a first order partial derivative
            % shall be calculated, is added to the supplied structure as a
            % field: param.X
            % X is from the set { 'R1', 'R2', 'D', 'B1', 'FlipAngle', 'Phase', 'tau', 'p', 's', 'S' }
            % the length of param.X specifies further selection
            % like specific tissue(s), RF pulse(s) (details see below)
            
            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            if ( isfield( param, 'R1' ) )
                
                cm.dR1 = param.R1;
                cm.len.dR1 = length( cm.dR1 );
                
                % in absence if spin coupling, the following definition (and equally for R2 and D) 
                % is a slight waste of allocated space (due to sparsity)
                % this is done in favor of a better readable (and maintainable) code
                
                cm.dm.dR1 = zeros( 3, cm.n_tissues, cm.alloc_n, 1, cm.len.dR1 ); 
                                                    
            end
            
            if ( isfield( param, 'R2' ) )
                
                cm.dR2 = param.R2;
                cm.len.dR2 = length( cm.dR2 );
                cm.dm.dR2 = zeros( 3, cm.n_tissues, cm.alloc_n, 1, cm.len.dR2 ); % cf. comment for R1
               
            end
            
            if ( isfield( param, 'D' ) )
                   
                cm.dD = param.D;
                cm.len.dD = length( cm.dD );
                cm.dm.dD = zeros( 3, cm.n_tissues, cm.alloc_n, cm.dim.D, cm.len.dD ); % cf. comment for R1
                cm.dEn_exp_dD = zeros( 3, 1, cm.alloc_n, cm.dim.D, cm.alloc_d ); 
                                
            end
            
            if ( isfield( param, 'B1' ) )
                
                cm.dB1 = 1;
                cm.len.dB1 = 1;
                cm.dm.dB1 = zeros( 3, cm.n_tissues, cm.alloc_n );
                
            end
            
            if ( isfield( param, 'FlipAngle' ) )
                
                cm.dFlipAngle = param.FlipAngle;
                cm.len.dFlipAngle = length( cm.dFlipAngle );
                cm.dm.dFlipAngle = zeros( 3, cm.n_tissues, cm.alloc_n, 1, cm.len.dFlipAngle );
                
            end
            
            if ( isfield( param, 'Phase' ) )
                    
                cm.dPhase = param.Phase;
                cm.len.dPhase = length( cm.dPhase );
                cm.dm.dPhase = zeros( 3, cm.n_tissues, cm.alloc_n, 1, cm.len.dPhase );
                    
            end
            
            if ( isfield( param, 'tau' ) )
                    
                cm.dtau = param.tau;
                cm.len.dtau = length( cm.dtau );
                cm.dm.dtau   = zeros( 3, cm.n_tissues, cm.alloc_n, 1, cm.len.dtau );
                cm.dEn_exp_dtau = zeros( 3, cm.n_tissues, cm.alloc_n, 1, cm.len.dtau );
                    
            end
                    
            if ( isfield( param, 'p' ) )

                cm.dp = param.p;
                cm.len.dp = length( cm.dp );
                cm.dm.dp   = zeros( 3, cm.n_tissues, cm.alloc_n, 3, cm.len.dp );
                cm.dEn_exp_dp = zeros( 3, cm.n_tissues, cm.alloc_n, 3, cm.len.dp, cm.alloc_d );
                
            end
                
            if ( isfield( param, 's' ) )

                if ( cm.diffusion == CoMoTk.No_Diff )
                    
                    error( 'Shape derivative without diffusion meaningless.' );
                    
                end

                cm.ds = param.s;
                cm.len.ds = length( cm.ds );
                cm.dm.ds   = zeros( 3, cm.n_tissues, cm.alloc_n, 3, cm.len.ds );
                cm.dEn_exp_ds = zeros( 1, cm.n_tissues, cm.alloc_n, 3, cm.len.ds );
                
            end
            
            if ( isfield( param, 'S' ) )
                
                if ( cm.diffusion == CoMoTk.No_Diff )
                    
                    error( 'Shape derivative without diffusion meaningless.' );
                    
                end
                
                cm.dS = param.S;
                cm.len.dS = length( cm.dS );
                cm.dm.dS   = zeros( 3, cm.n_tissues, cm.alloc_n, cm.dim.S, cm.len.dS );
                cm.dEn_exp_dS = zeros( 1, cm.n_tissues, 1, cm.dim.S, cm.len.dS );
                
            end
            
        end
        
        %% spin dynamics
        
        function RF ( cm, param )
            % Applies an instantaneous RF pulse
            %
            % IN
            %
            % param : structure with the info about the RF pulse,
            % specifically:
            %
            % param.FlipAngle : flip angle [rad] (mandatory parameter)
            %
            % param.Phase     : phase [rad] (mandatory parameter)
            % 
            % optional parameters to identify derivatives:
            %
            % param.handle_FlipAngle : handle to identify flip angle
            % 
            % param.handle_Phase     : handle to identify phase
            % 
            % param.MT_sat           : parameters for MT saturation pulse
            % (instead of instantaneous rotation) 
            % For proper use see example script and associated passage in the manuscript.
            
            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            % check for minimal information
            
            if ( ~isfield( param, 'FlipAngle' ) )
                
                error( 'No flip angle supplied.' );
                
            end
            
            if ( ~isfield( param, 'Phase' ) )
                
                error( 'No phase supplied.' );
                
            end
            
            % this is our actual flip angle
            
            ActualFlipAngle = cm.B1 * param.FlipAngle;
            
            % no action required, if no rotation is performed
            
            if ( mod( ActualFlipAngle, 2 * pi ) == 0 )
                
                return;
                
            end
            
            % check for derivative handles
            
            handle_FlipAngle = 0;
            handle_Phase = 0;
                
            if ( isfield( param, 'handle_FlipAngle' ) )
                
                handle_FlipAngle = param.handle_FlipAngle;
                
            end
            
            if ( isfield( param, 'handle_Phase' ) )
                
                handle_Phase = param.handle_Phase;
                
            end
            
            % calculate the rotation matrix
            
            c_al = cos( ActualFlipAngle );
            s_al = sin( ActualFlipAngle );
            
            ei_ph = cos( param.Phase ) + 1i * sin( param.Phase );
            
            RotMat = zeros( 3 );
            
            RotMat( 1, 1 ) = 0.5 * ( 1 + c_al );
            RotMat( 1, 2 ) = - 1i * s_al * ei_ph * CoMoTk.sqrt_0p5;
            RotMat( 1, 3 ) = 0.5 * ( 1 - c_al ) * ei_ph^2;
            RotMat( 2, 1 ) = - conj( RotMat( 1, 2 ) );
            RotMat( 2, 2 ) = c_al;
            RotMat( 2, 3 ) = - RotMat( 1, 2 );
            RotMat( 3, 1 ) = conj( RotMat( 1, 3 ) );
            RotMat( 3, 2 ) = - RotMat( 2, 1 );
            RotMat( 3, 3 ) = RotMat( 1, 1 );
            
            % check for magnetization transfer (MT)
            
            if ( isfield( param, 'MT_sat' ) )
            
                if ( length( param.MT_sat ) ~= cm.n_tissues )
                    
                    error( 'wrong size of saturation due to magnetization transfer' );
                    
                end
                
                % a value > 0 is treated as a bound pool
                
                b_MT = param.MT_sat > 0;
                
            else
                
                b_MT = false( cm.n_tissues, 1 );
                
            end
            
            n_MT = sum( b_MT );
            
            if ( n_MT > 0 )
                
                % Note: derivatives not (yet) implemented in presence of MT
                
                if ( cm.n_tissues > n_MT ) % should be always the case..
                
                    cm.m( :, ~b_MT, cm.b_n ) = ...
                        reshape( ...
                        RotMat * ...
                        reshape( ...
                        cm.m( :, ~b_MT, cm.b_n ), ...
                        [ 3, ( cm.n_tissues - n_MT ) * cm.n_conf ] ), ...
                        [ 3, ( cm.n_tissues - n_MT ), cm.n_conf ] );
                    
                end
 
                % MT saturation affects only the longitudinal component

                cm.m( 2, b_MT, cm.b_n ) = ...
                    reshape( param.MT_sat( b_MT ), [ 1, n_MT ] ) .* cm.m( 2, b_MT, cm.b_n );
                
                % since transverse magnetization is never generated, we
                % don't need to set it to zero explicitly

            else
                
                %% calculate derivatives, if necessary
                % This has to be done (at least for B1, FlipAngle and Phase) prior
                % to updating the configuration vector cm.m itself.
                % The rotation matrix does not depend on R1, R2, D, tau, p and s
                % which simplifes the associated recursions.
                
                % generally valid part of chain rule
                
                fn = fieldnames( cm.dim );
                
                for i = 1 : length( fn )
                    
                    X = fn{ i };
                    dX = [ 'd', X ];
                    
                    if ( cm.len.( dX ) > 0 )
                        
                        sz_tmp = [ 3, cm.n_tissues, cm.n_conf, cm.dim.( X ), cm.len.( dX ) ];
                        
                        % the following reshapes are required by the matrix multiplication
                        
                        cm.dm.( dX )( :, :, cm.b_n, :, : ) = ...
                            reshape( ...
                            RotMat * ...
                            reshape( ...
                            cm.dm.( dX )( :, :, cm.b_n, :, : ), ...
                            [ 3, prod( sz_tmp( 2 : end ) ) ] ), ...
                            sz_tmp );
                        
                    end
                    
                end
                
                % specific parts, when the derivative of the rotation matrix is nonzero
                
                % B1 derivative
                
                if ( cm.len.dB1 > 0 )
                    
                    dc_al = - s_al * param.FlipAngle;
                    ds_al = c_al * param.FlipAngle;
                    
                    % calculate the rotation matrix
                    
                    dRotMat_dB1 = zeros( 3 );
                    
                    dRotMat_dB1( 1, 1 ) = 0.5 * dc_al;
                    dRotMat_dB1( 1, 2 ) = - 1i * ds_al * ei_ph * CoMoTk.sqrt_0p5;
                    dRotMat_dB1( 1, 3 ) = - 0.5 * dc_al * ei_ph^2;
                    dRotMat_dB1( 2, 1 ) = - conj( dRotMat_dB1( 1, 2 ) );
                    dRotMat_dB1( 2, 2 ) = dc_al;
                    dRotMat_dB1( 2, 3 ) = - dRotMat_dB1( 1, 2 );
                    dRotMat_dB1( 3, 1 ) = conj( dRotMat_dB1( 1, 3 ) );
                    dRotMat_dB1( 3, 2 ) = - dRotMat_dB1( 2, 1 );
                    dRotMat_dB1( 3, 3 ) = dRotMat_dB1( 1, 1 );
                    
                    cm.dm.dB1( :, :, cm.b_n ) = cm.dm.dB1( :, :, cm.b_n ) + ...
                        reshape( ...
                        dRotMat_dB1 * ...
                        reshape( ...
                        cm.m( :, :, cm.b_n ), ...
                        [ 3, cm.n_tissues * cm.n_conf ] ), ...
                        [ 3, cm.n_tissues, cm.n_conf ] );
                    
                end
                
                % FlipAngle derivative(s)
                
                if ( cm.len.dFlipAngle > 0 )
                    
                    [ ~, i ] = find( cm.dFlipAngle == handle_FlipAngle );
                    
                    if ( ~isempty( i ) )
                        
                        dc_al = - s_al * cm.B1;
                        ds_al = c_al * cm.B1;
                        
                        % calculate the rotation matrix
                        
                        dRotMat_dFlipAngle = zeros( 3 );
                        
                        dRotMat_dFlipAngle( 1, 1 ) = 0.5 * dc_al;
                        dRotMat_dFlipAngle( 1, 2 ) = - 1i * ds_al * ei_ph * CoMoTk.sqrt_0p5;
                        dRotMat_dFlipAngle( 1, 3 ) = - 0.5 * dc_al * ei_ph^2;
                        dRotMat_dFlipAngle( 2, 1 ) = - conj( dRotMat_dFlipAngle( 1, 2 ) );
                        dRotMat_dFlipAngle( 2, 2 ) = dc_al;
                        dRotMat_dFlipAngle( 2, 3 ) = - dRotMat_dFlipAngle( 1, 2 );
                        dRotMat_dFlipAngle( 3, 1 ) = conj( dRotMat_dFlipAngle( 1, 3 ) );
                        dRotMat_dFlipAngle( 3, 2 ) = - dRotMat_dFlipAngle( 2, 1 );
                        dRotMat_dFlipAngle( 3, 3 ) = dRotMat_dFlipAngle( 1, 1 );
                        
                        cm.dm.dFlipAngle( :, :, cm.b_n, 1, i ) = cm.dm.dFlipAngle( :, :, cm.b_n, 1, i ) + ...
                            reshape( ...
                            dRotMat_dFlipAngle * ...
                            reshape( ...
                            cm.m( :, :, cm.b_n ), ...
                            [ 3, cm.n_tissues * cm.n_conf ] ), ...
                            [ 3, cm.n_tissues, cm.n_conf ] );
                        
                    end
                    
                end
                
                % Phase derivative(s)
                
                if ( cm.len.dPhase > 0 )
                    
                    [ ~, i ] = find( cm.dPhase == handle_Phase );
                    
                    if ( ~isempty( i ) )
                        
                        % calculate the rotation matrix
                        
                        dRotMat_dPhase = zeros( 3 );
                        
                        dRotMat_dPhase( 1, 2 ) = 1i * RotMat( 1, 2 );
                        dRotMat_dPhase( 1, 3 ) = 2 * 1i * RotMat( 1, 3 );
                        dRotMat_dPhase( 2, 1 ) = - conj( dRotMat_dPhase( 1, 2 ) );
                        dRotMat_dPhase( 2, 3 ) = - dRotMat_dPhase( 1, 2 );
                        dRotMat_dPhase( 3, 1 ) = conj( dRotMat_dPhase( 1, 3 ) );
                        dRotMat_dPhase( 3, 2 ) = - dRotMat_dPhase( 2, 1 );
                        
                        cm.dm.dPhase( :, :, cm.b_n, 1, i ) = cm.dm.dPhase( :, :, cm.b_n, 1, i ) + ...
                            reshape( ...
                            dRotMat_dPhase * ...
                            reshape( ...
                            cm.m( :, :, cm.b_n ), ...
                            [ 3, cm.n_tissues * cm.n_conf ] ), ...
                            [ 3, cm.n_tissues, cm.n_conf ] );
                        
                    end
                    
                end
            
                %% finally we update the configuration cm.m
                
                cm.m( :, :, cm.b_n ) = ...
                    reshape( ...
                    RotMat * ...
                    reshape( ...
                    cm.m( :, :, cm.b_n ), ...
                    [ 3, cm.n_tissues * cm.n_conf ] ), ...
                    [ 3, cm.n_tissues, cm.n_conf ] );

            end
                        
        end
        
        function b_is_ok = time ( cm, param )
            % Applies time interval
            %
            % IN
            %
            % param : structure with the info about the RF pulse,
            % specifically:
            %
            % param.lambda : unique index to identify the interval (mandatory parameter)
            % 
            % param.tau    : duration [ms] (mandatory in first call only)
            % 
            % param.p      : zero-order gradient moment [ rad / um ]
            % (only mandatory if D ~= 0)
            % 
            % param.s and param.S : information about gradient shape
            % (relevant for diffusion only, otherwise ignored)
            % 
            % param.v      : (constant) bulk motion velocity vector during
            % interval [ mm / s ] (optional)
            %
            % param.x      : location (vector) at start of interval [ um ]
            % (required in case of bulk motion)
            %
            % param.p1     : first order gradient moment [ rad * ms / um ]
            % (optional, only relevant in case of bulk motion)
            %
            % OUT
            %
            % b_is_ok : (bool) status, used for debugging
            
            % check for minimal information
            
            if ( ~isfield( param, 'lambda' ) )
                
                error( 'No ''lambda'' supplied.' );
                
            end

            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
                        
            b_is_ok = true;
            
            if ( cm.options.debug && ~cm.check_state( 'before pre_time' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            [ idx_time, b_n_new, b_up_free_time, b_do_free_time ] = cm.pre_time( param.lambda );
            
            if ( cm.options.debug && ~cm.check_state( 'after pre_time' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            % check optional arguments
            
            % is interval called for the first time?
            
            b_new_lambda = ~cm.b_tau( idx_time );

            % check for required gradient information in first call of time()
            
            if ( ...
                    b_new_lambda && ...
                    cm.diffusion ~= CoMoTk.No_Diff && ...
                    ~isfield( param, 'p' ) )
                
                error( 'Diffusion requires gradient moment in first call of ''time''.' );
                
            end            
            
            % duration
            
            if ( isfield( param, 'tau' ) )
                
                if ( b_new_lambda )
                    
                    cm.b_tau( idx_time ) = true;
                    cm.tau( idx_time ) = param.tau;
                    
                elseif ( cm.tau( idx_time ) ~= param.tau )
                    
                    error( 'Time interval tau must not change.' );
                    
                end
                
            elseif ( b_new_lambda )
                
                error( 'Missing duration in first call of ''time''.' );

            end
            
            % gradient moment
            
            if ( isfield( param, 'p' ) )
                
                if ( b_new_lambda )
                    
                    cm.b_p( idx_time ) = true;
                    cm.p( :, idx_time ) = param.p;
                    
                else
                    
                    if ( max( abs( cm.p( :, idx_time ) - param.p( : ) ) ) > 0 )
                        
                        error( 'Gradient moment p must not change.' );
                        
                    end
                    
                end
                
            end
                        
            if ( cm.diffusion ~= CoMoTk.No_Diff )
                
                % gradient shape (ignored, in absence of diffusion)
                
                if ( isfield( param, 's' ) )
                    
                    if ( ~b_new_lambda && ...
                            max( abs( cm.s( :, idx_time ) - param.s( : ) ) ) > 10 * eps )
                        
                        % Derivatives with respect to s make only sense for those intervals,
                        % for which the shape is not variable.
                        % This is checked by the following statement.
                        
                        [ ~, i_ ] = find( cm.ds == param.lambda );
                        
                        if ( ~isempty( i_ ) )
                            
                            error( 'Derivative of variable shape makes no sense.' );
                            
                        end
                        
                        % new shape invalidates previous damping factors
                        
                        cm.b_En( :, idx_time ) = false;
                        
                    end
                    
                    cm.s( :, idx_time ) = param.s( : );
                    
                elseif ( b_new_lambda )
                    
                    cm.s( :, idx_time ) = cm.p( :, idx_time );
                    
                end
                
                if ( isfield( param, 'S' ) )
                    
                    S_vec = zeros( cm.dim.S, 1 );
                    
                    if ( cm.diffusion == CoMoTk.Iso_Diff )
                        
                        S_vec( 1 ) = param.S;
                        
                    elseif ( cm.diffusion == CoMoTk.Tensor_Diff )
                        
                        S_vec( CoMoTk.Tensor_idx.xx ) = param.S( 1, 1 );
                        S_vec( CoMoTk.Tensor_idx.yy ) = param.S( 2, 2 );
                        S_vec( CoMoTk.Tensor_idx.zz ) = param.S( 3, 3 );
                        S_vec( CoMoTk.Tensor_idx.xy ) = param.S( 1, 2 );
                        S_vec( CoMoTk.Tensor_idx.xz ) = param.S( 1, 3 );
                        S_vec( CoMoTk.Tensor_idx.yz ) = param.S( 2, 3 );
                        
                    else
                        
                        error( 'This should not happen.' );
                        
                    end
                    
                    if ( ~b_new_lambda && ...
                            max( abs( cm.S( :, idx_time ) - S_vec( : ) ) ) > 10 * eps )
                        
                        % Derivatives with respect to s make only sense for those intervals,
                        % for which the shape is not variable.
                        % This is checked by the following statement.
                        
                        [ ~, i_ ] = find( cm.dS == param.lambda );
                        
                        if ( ~isempty( i_ ) )
                            
                            error( 'Derivative of variable shape makes no sense.' );
                            
                        end
                        
                        % new shape invalidates previous damping factors
                        
                        cm.b_En( :, idx_time ) = false;
                        
                    end

                    cm.S( :, idx_time ) = S_vec( : );
                    
                elseif ( b_new_lambda )
                    
                    S_vec = zeros( cm.dim.S, 1 );
                    
                    if ( cm.diffusion == CoMoTk.Iso_Diff )
                        
                        S_vec( 1 ) = ( cm.p( :, idx_time )' * cm.p( :, idx_time ) ) / 3;
                        
                    elseif ( cm.diffusion == CoMoTk.Tensor_Diff )
                        
                        S_tmp = ( cm.p( :, idx_time ) * cm.p( :, idx_time )' ) ./ 3;
                        
                        S_vec( CoMoTk.Tensor_idx.xx ) = S_tmp( 1, 1 );
                        S_vec( CoMoTk.Tensor_idx.yy ) = S_tmp( 2, 2 );
                        S_vec( CoMoTk.Tensor_idx.zz ) = S_tmp( 3, 3 );
                        S_vec( CoMoTk.Tensor_idx.xy ) = S_tmp( 1, 2 );
                        S_vec( CoMoTk.Tensor_idx.xz ) = S_tmp( 1, 3 );
                        S_vec( CoMoTk.Tensor_idx.yz ) = S_tmp( 2, 3 );
                        
                    else
                        
                        error( 'This should not happen.' );
                        
                    end
                    
                    cm.S( :, idx_time ) = S_vec( : );
                                        
                end
                
            end
            
            % update relaxation matrices
            
            cm.update_relaxation( idx_time );
            
            % chemical shift
            
            cs_tau = reshape( exp( - 1i .* cm.dom .* cm.tau( idx_time ) ), [ 1, cm.n_tissues ] );
                
            if ( isempty( cm.k ) )              % no magnetization exchange
                
                %% calculate recursive update of derivatives, if necessary
                % At least for R1, R2, D, tau, p and s, this has to be done prior
                % to updating the configuration vector cm.m itself.
                % The recursion coefficients do not depend on B1, FlipAngle and Phase
                % which simplifes the associated recursions.
                
                % generally valid part of chain rule
                
                fn = fieldnames( cm.dim );
                
                for i = 1 : length( fn )
                    
                    X = fn{ i };
                    dX = [ 'd', X ];
                    
                    if ( cm.len.( dX ) > 0 )
                        
                        if ( cm.diffusion == CoMoTk.No_Diff )
                            
                            cm.dm.( dX )( 1, :, cm.idx_up( cm.b_n, idx_time ), :, : ) = ...
                                cm.E( 2, :, idx_time ) .* ...
                                cm.dm.( dX )( 1, :, cm.b_n, :, : );
                            
                            cm.dm.( dX )( 2, :, cm.b_n, :, : ) = ...
                                cm.E( 1, :, idx_time ) .* ...
                                cm.dm.( dX )( 2, :, cm.b_n, :, : );
                            
                            cm.dm.( dX )( 3, :, cm.idx_do( cm.b_n, idx_time ), :, : ) = ...
                                cm.E( 2, :, idx_time ) .* ...
                                cm.dm.( dX )( 3, :, cm.b_n, :, : );
                            
                        else
                            
                            cm.dm.( dX )( 1, :, cm.idx_up( cm.b_n, idx_time ), :, : ) = ...
                                cm.En( 1, :, cm.b_n, idx_time ) .* ...
                                cm.dm.( dX )( 1, :, cm.b_n, :, : );
                            
                            cm.dm.( dX )( 2, :, cm.b_n, :, : ) = ...
                                cm.En( 2, :, cm.b_n, idx_time ) .* ...
                                cm.dm.( dX )( 2, :, cm.b_n, :, : );
                            
                            cm.dm.( dX )( 3, :, cm.idx_do( cm.b_n, idx_time ), :, : ) = ...
                                cm.En( 3, :, cm.b_n, idx_time ) .* ...
                                cm.dm.( dX )( 3, :, cm.b_n, :, : );
                            
                        end
                        
                        % cleanup of free transverse configurations
                        
                        cm.dm.( dX )( 1, :, b_do_free_time, :, : ) = 0;
                        cm.dm.( dX )( 3, :, b_up_free_time, :, : ) = 0;
                        
                    end
                    
                end
                
                % specific parts, when the derivative of the damping matrix and/or the recovery term is nonzero
                
                % R1 derivative(s)
                
                if ( cm.len.dR1 > 0 )
                    
                    for i = 1 : cm.len.dR1
                        
                        if ( cm.diffusion == CoMoTk.No_Diff )
                            
                            cm.dm.dR1( 2, cm.dR1( i ), cm.b_n, 1, i ) = ...
                                cm.dm.dR1( 2, cm.dR1( i ), cm.b_n, 1, i ) - ...
                                cm.tau( idx_time ) .* cm.E( 1, cm.dR1( i ), idx_time ) .* cm.m( 2, cm.dR1( i ), cm.b_n );
                            
                        else
                            
                            cm.dm.dR1( 2, cm.dR1( i ), cm.b_n, 1, i ) = ...
                                cm.dm.dR1( 2, cm.dR1( i ), cm.b_n, 1, i ) - ...
                                cm.tau( idx_time ) .* cm.En( 2, cm.dR1( i ), cm.b_n, idx_time ) .* cm.m( 2, cm.dR1( i ), cm.b_n ); ...
                                
                        end
                        
                        cm.dm.dR1( 2, cm.dR1( i ), cm.null_idx, 1, i ) = ...
                            cm.dm.dR1( 2, cm.dR1( i ), cm.null_idx, 1, i ) + ...
                            cm.mu( cm.dR1( i ) ) .* cm.tau( idx_time ) .* cm.E( 1, cm.dR1( i ), idx_time );
                        
                    end
                    
                end
                
                % R2 derivative(s)
                
                if ( cm.len.dR2 > 0 )
                    
                    for i = 1 : cm.len.dR2
                        
                        if ( cm.diffusion == CoMoTk.No_Diff )
                            
                            cm.dm.dR2( 1, cm.dR2( i ), cm.idx_up( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dR2( 1, cm.dR2( i ), cm.idx_up( cm.b_n, idx_time ), 1, i ) - ...
                                cm.tau( idx_time ) .* cs_tau( cm.dR2( i ) ) .* cm.E( 2, cm.dR2( i ), idx_time ) .* cm.m( 1, cm.dR2( i ), cm.b_n );
                            
                            cm.dm.dR2( 3, cm.dR2( i ), cm.idx_do( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dR2( 3, cm.dR2( i ), cm.idx_do( cm.b_n, idx_time ), 1, i ) - ...
                                cm.tau( idx_time ) .* conj( cs_tau( cm.dR2( i ) ) ) .* cm.E( 2, cm.dR2( i ), idx_time ) .* cm.m( 3, cm.dR2( i ), cm.b_n );
                            
                        else
                            
                            cm.dm.dR2( 1, cm.dR2( i ), cm.idx_up( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dR2( 1, cm.dR2( i ), cm.idx_up( cm.b_n, idx_time ), 1, i ) - ...
                                cm.tau( idx_time ) .* cs_tau( cm.dR2( i ) ) .* cm.En( 1, cm.dR2( i ), cm.b_n, idx_time ) .* cm.m( 1, cm.dR2( i ), cm.b_n );
                            
                            cm.dm.dR2( 3, cm.dR2( i ), cm.idx_do( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dR2( 3, cm.dR2( i ), cm.idx_do( cm.b_n, idx_time ), 1, i ) - ...
                                cm.tau( idx_time ) .* conj( cs_tau( cm.dR2( i ) ) ) .* cm.En( 3, cm.dR2( i ), cm.b_n, idx_time ) .* cm.m( 3, cm.dR2( i ), cm.b_n );
                            
                        end
                        
                    end
                    
                end
                
                % D derivative(s)
                
                if ( cm.len.dD > 0 )
                    
                    for i = 1 : cm.len.dD
                        
                        cm.dm.dD( 1, cm.dD( i ), cm.idx_up( cm.b_n, idx_time ), :, i ) = ...
                            cm.dm.dD( 1, cm.dD( i ), cm.idx_up( cm.b_n, idx_time ), :, i ) + ...
                            cs_tau( cm.dD( i ) ) .* ...
                            cm.dEn_exp_dD( 1, 1, cm.b_n, :, idx_time ) .* ...
                            cm.En( 1, cm.dD( i ), cm.b_n, idx_time ) .* ...
                            cm.m( 1, cm.dD( i ), cm.b_n );
                        
                        cm.dm.dD( 2, cm.dD( i ), cm.b_n, :, i ) = ...
                            cm.dm.dD( 2, cm.dD( i ), cm.b_n, :, i ) + ...
                            cm.dEn_exp_dD( 2, 1, cm.b_n, :, idx_time ) .* ...
                            cm.En( 2, cm.dD( i ), cm.b_n, idx_time ) .* ...
                            cm.m( 2, cm.dD( i ), cm.b_n );
                        
                        cm.dm.dD( 3, cm.dD( i ), cm.idx_do( cm.b_n, idx_time ), :, i ) = ...
                            cm.dm.dD( 3, cm.dD( i ), cm.idx_do( cm.b_n, idx_time ), :, i ) + ...
                            conj( cs_tau( cm.dD( i ) ) ) .* ...
                            cm.dEn_exp_dD( 3, 1, cm.b_n, :, idx_time ) .* ...
                            cm.En( 3, cm.dD( i ), cm.b_n, idx_time ) .* ...
                            cm.m( 3, cm.dD( i ), cm.b_n );
                        
                    end
                    
                end
                
                % tau derivative(s)
                
                if ( cm.len.dtau > 0 )
                    
                    [ ~, i ] = find( cm.dtau == param.lambda );
                    
                    if ( ~isempty( i ) )
                        
                        if ( cm.diffusion == CoMoTk.No_Diff )
                            
                            cm.dm.dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), 1, i ) - ...
                                cm.R2 .* cs_tau .* cm.E( 2, :, idx_time ) .* cm.m( 1, :, cm.b_n );
                            
                            cm.dm.dtau( 2, :, cm.b_n, 1, i ) = ...
                                cm.dm.dtau( 2, :, cm.b_n, 1, i ) - ...
                                cm.R1 .* cm.E( 1, :, idx_time ) .* cm.m( 2, :, cm.b_n );
                            
                            cm.dm.dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), 1, i ) - ...
                                cm.R2 .* conj( cs_tau ) .* cm.E( 2, :, idx_time ) .* cm.m( 3, :, cm.b_n );
                            
                        else
                            
                            cm.dm.dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dtau( 1, :, cm.idx_up( cm.b_n, idx_time ), 1, i ) + ...
                                ( cm.dEn_exp_dtau( 1, :, cm.b_n, 1, i ) - cm.R2 ) .* ...
                                cs_tau .* cm.En( 1, :, cm.b_n, idx_time ) .* cm.m( 1, :, cm.b_n );
                            
                            cm.dm.dtau( 2, :, cm.b_n, 1, i ) = ...
                                cm.dm.dtau( 2, :, cm.b_n, 1, i ) + ...
                                ( cm.dEn_exp_dtau( 2, :, cm.b_n, 1, i ) - cm.R1 ).* ...
                                cm.En( 2, :, cm.b_n, idx_time ) .* cm.m( 2, :, cm.b_n );
                            
                            cm.dm.dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), 1, i ) = ...
                                cm.dm.dtau( 3, :, cm.idx_do( cm.b_n, idx_time ), 1, i ) + ...
                                ( cm.dEn_exp_dtau( 3, :, cm.b_n, 1, i ) - cm.R2 ) .* ...
                                conj( cs_tau ) .* cm.En( 3, :, cm.b_n, idx_time ) .* cm.m( 3, :, cm.b_n );
                            
                        end
                        
                        cm.dm.dtau( 2, :, cm.null_idx, 1, i ) = cm.dm.dtau( 2, :, cm.null_idx, 1, i ) + ...
                            cm.mu .* cm.R1 .* cm.E( 1, :, idx_time );
                        
                    end
                    
                end
                
                % p derivative(s)
                
                if ( cm.len.dp > 0 )
                    
                    cm.dm.dp( 1, :, cm.idx_up( cm.b_n, idx_time ), :, : ) = ...
                        cm.dm.dp( 1, :, cm.idx_up( cm.b_n, idx_time ), :, : ) + ...
                        cs_tau .* ...
                        cm.dEn_exp_dp( 1, :, cm.b_n, :, :, idx_time ) .* ...
                        cm.En( 1, :, cm.b_n, idx_time ) .* ...
                        cm.m( 1, :, cm.b_n );
                    
                    cm.dm.dp( 2, :, cm.b_n, :, : ) = ...
                        cm.dm.dp( 2, :, cm.b_n, :, : ) + ...
                        cm.dEn_exp_dp( 2, :, cm.b_n, :, :, idx_time ) .* ...
                        cm.En( 2, :, cm.b_n, idx_time ) .* ...
                        cm.m( 2, :, cm.b_n );
                    
                    cm.dm.dp( 3, :, cm.idx_do( cm.b_n, idx_time ), :, : ) = ...
                        cm.dm.dp( 3, :, cm.idx_do( cm.b_n, idx_time ), :, : ) + ...
                        conj( cs_tau ) .* cm.dEn_exp_dp( 3, :, cm.b_n, :, :, idx_time ) .* ...
                        cm.En( 3, :, cm.b_n, idx_time ) .* ...
                        cm.m( 3, :, cm.b_n );
                    
                end
                
                % s derivative(s)
                
                if ( cm.len.ds > 0 )
                    
                    [ ~, i ] = find( cm.ds == param.lambda );
                    
                    if ( ~isempty( i ) )
                        
                        cm.dm.ds( 1, :, cm.idx_up( cm.b_n, idx_time ), :, i ) = ...
                            cm.dm.ds( 1, :, cm.idx_up( cm.b_n, idx_time ), :, i ) + ...
                            cs_tau .* cm.dEn_exp_ds( 1, :, cm.b_n, :, i ) .* cm.En( 1, :, cm.b_n, idx_time ) .* cm.m( 1, :, cm.b_n );
                        
                        cm.dm.ds( 3, :, cm.idx_do( cm.b_n, idx_time ), :, i ) = ...
                            cm.dm.ds( 3, :, cm.idx_do( cm.b_n, idx_time ), :, i ) - ...
                            conj( cs_tau ) .* cm.dEn_exp_ds( 1, :, cm.b_n, :, i ) .* cm.En( 3, :, cm.b_n, idx_time ) .* cm.m( 3, :, cm.b_n );
                        
                    end
                    
                end
                
                % S derivative(s)
                
                if ( cm.len.dS > 0 )
                    
                    [ ~, i ] = find( cm.dS == param.lambda );
                    
                    if ( ~isempty( i ) )
                        
                        cm.dm.dS( 1, :, cm.idx_up( cm.b_n, idx_time ), :, i ) = ...
                            cm.dm.dS( 1, :, cm.idx_up( cm.b_n, idx_time ), :, i ) + ...
                            cs_tau .* cm.dEn_exp_dS( 1, :, 1, :, i ) .* cm.En( 1, :, cm.b_n, idx_time ) .* cm.m( 1, :, cm.b_n );
                        
                        cm.dm.dS( 3, :, cm.idx_do( cm.b_n, idx_time ), :, i ) = ...
                            cm.dm.dS( 3, :, cm.idx_do( cm.b_n, idx_time ), :, i ) + ...
                            conj( cs_tau ) .* cm.dEn_exp_dS( 1, :, 1, :, i ) .* cm.En( 3, :, cm.b_n, idx_time ) .* cm.m( 3, :, cm.b_n );
                        
                    end
                    
                end
                
                %% finally we update the configuration cm.m
                
                % relaxation
                
                if ( cm.diffusion == CoMoTk.No_Diff )
                    
                    cm.m( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                        cs_tau .* cm.E( 2, :, idx_time ) .* cm.m( 1, :, cm.b_n );
                    
                    cm.m( 2, :, cm.b_n ) = ...
                        cm.E( 1, :, idx_time ) .* cm.m( 2, :, cm.b_n );
                    
                    cm.m( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                        conj( cs_tau ) .* cm.E( 2, :, idx_time ) .* cm.m( 3, :, cm.b_n );
                    
                else
                    
                    cm.m( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                        cs_tau .* cm.En( 1, :, cm.b_n, idx_time ) .* cm.m( 1, :, cm.b_n );
                    
                    cm.m( 2, :, cm.b_n ) = ...
                        cm.En( 2, :, cm.b_n, idx_time ) .* cm.m( 2, :, cm.b_n );
                    
                    cm.m( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                        conj( cs_tau ) .* cm.En( 3, :, cm.b_n, idx_time ) .* cm.m( 3, :, cm.b_n );
                    
                end
                
                cm.m( 1, :, b_do_free_time ) = 0;
                cm.m( 3, :, b_up_free_time ) = 0;
                
                % repolarization (here, the relative proton densities come into play)
                
                cm.m( 2, :, cm.null_idx ) = cm.m( 2, :, cm.null_idx ) + cm.mu .* ( 1 - cm.E( 1, :, idx_time ) );
                
            else                              % with magnetization exchange

                %% magnetization exchange
                
                b_inf_R2 = isinf( cm.R2 );
                
                if ( any( b_inf_R2 ) ) 
                    
                    % In case of MT, there is no coupling of transverse
                    % magnetization to bound pool(s).
                    %
                    % For the (rare?) case than there is more than one
                    % mobile compartment, we allow for coupling between
                    % them, though. To this end, we define a transition matrix
                    % for the free compartments only and modify the sum
                    % rule accordingly (otherwise the diagonal elements
                    % would leak into the bound pool, which is usually
                    % assumed as uncoupled with respect to the transverse
                    % magnetization).
                    %
                    % NOTE: 
                    % Use with care, since I am not sure, whether this
                    % approach is sound.
                    
                    k_ = cm.k;
                    
                    % eliminate coupling to bound pool(s)
                    
                    k_( :, b_inf_R2 ) = 0;
                    k_( b_inf_R2, : ) = 0;
                    
                    % modify sum rule for transverse components
                    
                    for j = 1 : cm.n_tissues
                        
                        k_( j, j ) = 0;
                        
                        for i = 1 : cm.n_tissues
                            
                            if ( i ~= j )
                                
                                k_( j, j ) = k_( j, j ) - k_( i, j );
                                
                            end
                            
                        end
                        
                    end

                    Y_1 = k_ - diag( 1i .* cm.dom( : ) + cm.R2( : ) );
                    
                    Y_1( :, b_inf_R2 ) = 0;
                    Y_1( b_inf_R2, : ) = 0;
                    
                else
                    
                    Y_1 = cm.k - diag( 1i .* cm.dom( : ) + cm.R2( : ) );

                end
                
                Y_0 = cm.k - diag( cm.R1( : ) );

                exp_Y_tau_1 = expm( Y_1 .* cm.tau( idx_time ) );
                exp_Y_tau_0 = expm( Y_0 .* cm.tau( idx_time ) );

                y = cm.R1( : ) .* cm.mu( : );
                
                % relaxation
                
                cm.m( 1, :, cm.idx_up( cm.b_n, idx_time ) ) = ...
                    exp_Y_tau_1 * ...
                    reshape( ...
                    cm.m( 1, :, cm.b_n ), ...
                    [ cm.n_tissues, cm.n_conf ] );
                
                cm.m( 2, :, cm.b_n ) = ...7
                    exp_Y_tau_0 * ...
                    reshape( ...
                    cm.m( 2, :, cm.b_n ), ...
                    [ cm.n_tissues, cm.n_conf ] );
                
                cm.m( 3, :, cm.idx_do( cm.b_n, idx_time ) ) = ...
                    conj( exp_Y_tau_1 ) * ...
                    reshape( ...
                    cm.m( 3, :, cm.b_n ), ...
                    [ cm.n_tissues, cm.n_conf ] );
                
                cm.m( 1, :, b_do_free_time ) = 0;
                cm.m( 3, :, b_up_free_time ) = 0;
                
                % repolarization (here, the relative proton densities come into play)
                
                cm.m( 2, :, cm.null_idx ) = cm.m( 2, :, cm.null_idx ) - ...
                    reshape( ...
                    Y_0 \ ( ( eye( cm.n_tissues ) - exp_Y_tau_0 ) * y ), ...
                    [ 1, cm.n_tissues ] );
                
            end
            
            % perform final bookkeeping                        
            
            cm.post_time( b_n_new );
            
            if ( cm.options.debug && ~cm.check_state( 'after post_time' ) )
                
                b_is_ok = false;
                return;
                
            end
           
            % remove negligible configurations
            
            if( ~cm.meltdown )
                
                b_is_ok = false;
                return;
                
            end
            
            if ( cm.options.debug && ~cm.check_state( 'after meltdown' ) )
                
                b_is_ok = false;
                return;
                
            end
            
            % bulk motion due to constant velocity during the interval
            
            if ( isfield( param, 'v' ) )
                
                % "v" is the actual velocity vector, units are [ mm / s ]
                % (it does not matter, whether it is a row or column vector)

                if ( ~isfield( param, 'x' ) )
                    
                    error( 'The location should be supplied in case of bulk motion.' );
                    
                end                
                
                if ( isfield( param, 'p1' ) )
                    
                    % If the first order gradient is supplied, we use it.
                    % Unit of p1: [ rad * ms / um ]

                    dth = - param.p1( : )' * param.v( : );

                else
                    
                    % 0therwise, we assume a constant gradient during the
                    % interval. In that case, the angle due to dephasing is just
                    % - (tau/2) * (p * v)
                    
                    % NOTE: nonzero velocities are currently not taken into account, when calculating derivatives!
                    
                    dth = - ( 0.5 * cm.tau( idx_time ) ) * ( cm.p( :, idx_time )' * param.v( : ) );
                    
                end
                
                dth = dth - cm.p( :, idx_time )' * param.x( : );
                
                ei_dth = exp( 1i * dth );
                
                cm.m( 1, cm.b_n ) = ei_dth .* cm.m( 1, cm.b_n );
                cm.m( 3, cm.b_n ) = conj( ei_dth ) .* cm.m( 3, cm.b_n );
                
            end
            
            % update configuration dependent time and gradient
            
            cm.update_tau_n;
            cm.update_p_n;
            
        end
        
        function spoiler ( cm )
            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            % set all transverse magnetization to zero (instantaneously)
           
            cm.m( [ 1, 3 ], :, cm.b_n ) = 0;
            
            % remove negligible configurations
            
            cm.meltdown;
            
        end
        
        %% get results
        
        function b_n = find ( cm, lambda, n )
            % search for specific configurations
            %
            % IN
            %
            % lambda : array of CM dimension indices (== mandatory
            % parameter in method 'time')
            % n      : array of configuration orders
            %
            % OUT
            %
            % b_n : location of stored configurations, which match ANY of
            % the supplied condition pairs lambda(i), n(i)
            
            if ( length( lambda ) ~= length( n ) )
                
                error( 'lambda and n must have same length.' );
                
            end
            
            % check, which of the supplied lambda are actually alive (i.e. stored)
            
            b_lambda_2d = cm.b_lambda & ( lambda( : ) == cm.lambda );  % size( b_lambda_2d ) == [ length( lambda ), cm.alloc_d ]
            
            % non-existing lambda must not have non-zero configuration orders
            
            b_lambda_found = any( b_lambda_2d, 2 );            % size( b_lambda_found ) == [ length( lambda ), 1 ]
            
            if ( any( n( ~b_lambda_found ) ) )
                
                b_n = [];
                return;
                
            end

            % n == 0 is always true for not (yet) occupied time index lambda
            
            if ( sum( ~b_lambda_found ) > 0 )

                b_n = cm.b_n;
                return;
                
            end
            
            n_lambda_found = sum( b_lambda_found );
            
            % for the rest, determine the indices in storage
            
            [ idx_lambda, ~ ] = find( b_lambda_2d' );
            
            % we can now determine b_n
            
            b_n = cm.b_n & any( cm.n( :, idx_lambda ) == reshape( n( b_lambda_found ), [ 1, n_lambda_found ] ), 2 );
            
        end
        
        function res = sum ( cm, param )            
            % return the weighted CM sum and calculated derivatives
            %
            % IN
            %
            % structure param (if empty, complete sum over all stored
            % configurations is returned)
            %
            % param.b_n : over which configurations shall the sum be taken?
            % (optional)
            % (if field is missing, the sum is over all stored
            % configurations)
            %
            % param.omega : bulk off-resonance [ rad / ms ] (optional)
            % (default == 0, if field is missing)
            %
            % param.x : position vector [ um ] (optional)
            % (== zero vector, if field is missing)
            %
            % param.w_n : external weighting factors (optional, if supplied, 
            % size needs to match the number of addends in the CM sum)
            %
            % decay due to inhomogeneous broadening is taken into account,
            % if cm.inhomogeneous_decay is nonempty
            % (more detailed information about the format of this function
            % can be found above at the associated property
            %
            % OUT
            %
            % res.xy : transverse component (conventional transverse
            % magnetization)
            %
            % res.z  : longitudinal component
            %
            % res.dxy.dX : derivative of transverse component with respect
            % to X
            %
            % res.dx.dX : same for longitudinal component
            
            if ( ~isempty( param ) && isfield( param, 'b_n' ) )
                
                b_n_ = param.b_n;
                
            else
                
                b_n_ = cm.b_n;
                
            end
            
            if ( ~isempty( param ) && isfield( param, 'omega' ) )
                
                omega_ = param.omega;
                
            else
                
                omega_ = 0;
                
            end
            
            n_c = sum( b_n_ );
            
            if ( n_c == 0 )
                
                res = [];
                return;
                
            end
            
            % off-resonance and chemical shift (if applicable )
            
            w_n_ = exp( - 1i .* omega_ .* reshape( cm.tau_n( b_n_ ), [ 1, n_c ] ) );
            
            % *= gradient effects (optional)
            
            if ( ~isempty( param ) && isfield( param, 'x' ) )
                
                w_n_ = w_n_ .* exp ( - 1i .* sum( cm.p_n( :, b_n_ ) .* param.x( : ), 1 ) );
                                                                                          
            end

            % *= effects due to inhomogeneous broadening (if applicable, most commonly due to R2p)
            
            if ( ~isempty( cm.inhomogeneous_decay ) )
                
                w_n_ = w_n_ .* cm.inhomogeneous_decay( cm.tau_n( b_n_ ), param );
                
            end
            
            % *= explicit weights (optional)
            
            if ( ~isempty( param ) && isfield( param, 'w_n' ) )
                
                w_n_( : ) = w_n_( : ) .* param.w_n( : );
                                                                                          
            end

            % add a singleton dimension for coordinates
            
            sz_w_n = size( w_n_ );
            w_n_ = reshape( w_n_, [ 1, sz_w_n ] );
            
            % calculate the isochromat
            
            res.xy = CoMoTk.sqrt_2 .* sum( sum( w_n_ .* cm.m( 1, :, b_n_ ), 2 ), 3 );
            res.z = sum( sum( w_n_ .* cm.m( 2, :, b_n_ ), 2 ), 3 );  % should be real...
            
            % and the calculated derivatives
            
            fn = fieldnames( cm.dim );
            
            for i = 1 : length( fn )
                
                X = fn{ i };
                dX = [ 'd', X ];
                
                if ( cm.len.( dX ) > 0 )
                    
                    res.dxy.( dX ) = reshape( ...
                        CoMoTk.sqrt_2 .* sum( sum( ...
                        w_n_ .* cm.dm.( dX )( 1, :, b_n_, :, : ) ...
                        , 2 ), 3 ), ...
                        [ cm.dim.( X ), cm.len.( dX ) ] );
                    
                    res.dz.( dX ) = reshape( ...
                        sum( sum( ...
                        w_n_ .* cm.dm.( dX )( 2, :, b_n_, :, : ) ...
                        , 2 ), 3 ), ...
                        [ cm.dim.( X ), cm.len.( dX ) ] );
                    
                end
                
            end
                        
        end
        
        %% Lorentz function
        
        function res = lorentz_decay ( cm, tau, param )
            % calculates Lorentz decay exp( - tau * R2p )
            %
            % param.R2p : relaxation rate R2^\prime [ 1 / ms ] (mandatory parameter)
           
            if ( isfield( param, 'R2p' ) )
                
                R2p = param.R2p;
                
                if ( length( R2p ) ~= cm.n_tissues && length( R2p ) ~= 1 )
                    
                    error( 'length( R2p ) must be equal to the number of tissues or one.' );
                    
                end
                
                if ( length( R2p ) == cm.n_tissues )
                    
                    R2p = reshape( R2p, [ cm.n_tissues, 1 ] );
                    
                else
                    
                    R2p = repmat( R2p, [ cm.n_tissues, 1 ] );
                    
                end
                                
                res = exp( - R2p .* reshape( abs( tau ), [ 1, length( tau ) ] ) );
                               
            else
                
                error( 'R2p not specified in lorentz_decay' );
                
            end
            
        end
        
        function res = spherical_decay ( cm, tau, param )
            % decay due to spherical frequency distribution
            %
            % param.bw_sphere : max. off-resonance frequency magnitude 
            % [ rad / ms ] (mandatory parameter)
            
            if ( isfield( param, 'bw_sphere' ) )
                
                sigma = 0.5 .* param.bw_sphere;
                
                if ( length( sigma ) ~= cm.n_tissues && length( sigma ) ~= 1 )
                    
                    error( 'length( bw_sphere ) must be equal to the number of tissues or one.' );
                    
                end
                
                if ( length( sigma ) == cm.n_tissues )
                    
                    sigma = reshape( sigma, [ cm.n_tissues, 1 ] );
                    
                else
                    
                    sigma = repmat( sigma, [ cm.n_tissues, 1 ] );
                    
                end
                                    
                st_ = sigma .* reshape( tau, [ 1, length( tau ) ] );
                    
                res = besselj( 1, 2 .* st_ ) ./ ( st_ );

                res( st_ == 0 ) = 1;
                    
            else
                
                error( 'bw_sphere not specified in spherical_decay' );
                
            end
            
        end
        
        %% various checks for debugging purposes
        
        function b_is_ok = check_state ( cm, context )
            % used for debugging purposes
            
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
                
                cm.log.b_lambda{ cm.log.count } = cm.b_lambda;
                cm.log.lambda{ cm.log.count } = cm.lambda;
                
                cm.log.b_up_free{ cm.log.count } = cm.b_up_free;
                cm.log.b_do_free{ cm.log.count } = cm.b_do_free;
                
                cm.log.idx_up{ cm.log.count } = cm.idx_up;
                cm.log.idx_do{ cm.log.count } = cm.idx_do;
                
            end
            
            if ( cm.options.debug )
                
                b_up_occ = cm.b_n & ~cm.b_up_free & cm.b_lambda;
                b_do_occ = cm.b_n & ~cm.b_do_free & cm.b_lambda;
                
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
                
                if ( sum( b_ok ) ~= length( idx_up_diff ) )
                    
                    b_is_ok = false;
                    
                end
                
                idx_do_diff = cm.idx_up( sub2ind( [ cm.alloc_n, cm.alloc_d ], cm.idx_do( b_do_occ ), c_do ) ) - idx_do_occ;
                b_ok = idx_do_diff == 0;
                
                if ( cm.options.debug )
                    
                    cm.log.idx_do_diff{ cm.log.count } = idx_do_diff;
                    cm.log.b_idx_do_diff{ cm.log.count } = b_ok;
                    
                end
                
                if ( sum( b_ok ) ~= length( idx_do_diff ) )
                    
                    b_is_ok = false;
                    
                end
                
                n_diff = cm.n( sub2ind( [ cm.alloc_n, cm.alloc_d ], cm.idx_up( b_up_occ ), c_up ) ) - cm.n( b_up_occ );
                b_ok = n_diff == 1;
                
                if ( cm.options.debug )
                    
                    cm.log.n_diff{ cm.log.count } = n_diff;
                    cm.log.b_n_diff{ cm.log.count } = b_ok;
                    
                end
                
                if ( sum( b_ok ) ~= length( n_diff ) )
                    
                    b_is_ok = false;
                    
                end
                
            end
            
        end        
        
    end
    
    methods ( Access = private )
        % several helper routines

        function check ( cm )
            % Check, whether all mandatory tissue parameters have been set.
            
            if ( isempty( cm.R1 ) || ...
                    isempty( cm.R2 ) || ...
                    isempty( cm.D ) )
                
                error( 'Incomplete specification of mandatory tissue parameters.' );
                
            end
            
            % check proton densities
            
            if ( size( cm.mu, 2 ) ~= cm.n_tissues )
                
                error( 'Proton density not set properly.' );
                
            end
            
            % check correct long time limit of magnetization transfer/exchange
            % (particle conservation is properly incorporated by set.k) 
            
            if ( ~isempty( cm.k ) )
               
                % cm.mu must lie in the null space of cm.k:
                
                k_mu = cm.k * cm.mu( : );
                
                if ( max( abs( k_mu ) > 10 * eps ) )
                    
                    % Maybe there is an error, maybe it is just numerical
                    % inaccuracy. Provide an informal message to the user and let
                    % the program continue.
                    
                    fprintf( 1, 'Check, whether cm.mu lies in null space of cm.k\n' );
                    fprintf( 1, 'Found result for cm.k * cm.mu =\n' );
                    disp( k_mu );
                    
                end
                
            end
            
        end
        
        function prep ( cm )
            % Initializes magnetization and other stuff
            
            % check, whether required parameters have been set properly
            
            cm.check;
            
            % no diffusion AND magnetization exchange
            
            if ( cm.diffusion ~= CoMoTk.No_Diff && ~isempty( cm.k ) )
                
                error( 'Diffusion AND magnetization exchange not supported.')
                
            end
            
            % take initial settings from options
            
            cm.alloc_d = cm.options.alloc_d;
            cm.alloc_n = cm.options.alloc_n;
            cm.epsilon = cm.options.epsilon;
            cm.rapid_meltdown = cm.options.rapid_meltdown;
            cm.verbose = cm.options.verbose;
            cm.debug = cm.options.debug;
            
            % generate initial configuration vector (thermal equlibrium)
            
            cm.m = zeros( 3, cm.n_tissues, cm.alloc_n );
            cm.m( 2, :, cm.null_idx ) = cm.mu;
            
            % allocate space
            
            cm.b_lambda = false( 1, cm.alloc_d );
            cm.b_n = false( cm.alloc_n, 1 );
            
            cm.lambda = zeros( 1, cm.alloc_d );
            cm.n = zeros( cm.alloc_n, cm.alloc_d );
            
            cm.b_up_free = true( cm.alloc_n, cm.alloc_d );
            cm.b_do_free = true( cm.alloc_n, cm.alloc_d );
            
            cm.idx_up = zeros( cm.alloc_n, cm.alloc_d );
            cm.idx_do = zeros( cm.alloc_n, cm.alloc_d );
            
            cm.idx_up_new = zeros( cm.alloc_n, 1 );
            cm.idx_do_new = zeros( cm.alloc_n, 1 );
            
            % pure relaxation
            
            cm.b_E = false( cm.alloc_d, 1 );
            cm.E = zeros( 2, cm.n_tissues, cm.alloc_d );

            % diffusion
            
            if ( cm.diffusion ~= CoMoTk.No_Diff )
                
                cm.b_En = false( cm.alloc_n, cm.alloc_d );
                cm.En = zeros( 3, cm.n_tissues, cm.alloc_n, cm.alloc_d );
                cm.En_exp = zeros( 3, cm.n_tissues, cm.alloc_n, cm.alloc_d );
                
                if ( cm.diffusion == CoMoTk.Iso_Diff )
                    
                    cm.dim.D = 1;                    
                    cm.dim.S = 1;                    
                    
                elseif ( cm.diffusion == CoMoTk.Tensor_Diff )
                    
                    cm.dim.D = 6;                    
                    cm.dim.S = 6;                    
                    
                else
                    
                    error( 'This should not happen.' );
                    
                end

                cm.s = zeros( 3, cm.alloc_d );
                cm.S = zeros( cm.dim.S, cm.alloc_d );

            end
            
            % only the zeroth configuration is occupied
            
            cm.b_n( cm.null_idx ) = true;
            cm.n_conf = 1;
            
            % initialize times and gradients
            
            cm.b_tau = false( cm.alloc_d, 1 );
            cm.tau = zeros( cm.alloc_d, 1 );
            
            cm.b_p = false( cm.alloc_d, 1 );
            cm.p = zeros( 3, cm.alloc_d );
            
            cm.b_tau_n = false( cm.alloc_n, 1 );
            cm.tau_n = zeros( cm.alloc_n, 1 );

            cm.b_p_n = false( cm.alloc_n, 1 );
            cm.p_n = zeros( 3, cm.alloc_n );
                        
            % everything is prepared now
            
            cm.b_prep = true;
            
        end        
        
        function update_tau_n ( cm )
            % update cm.tau_n and cm.b_tau_n
            
            b_todo = cm.b_n & ~cm.b_tau_n;
            
            n_todo = sum( b_todo );
            
            if ( n_todo > 0 )
                
                cm.tau_n( b_todo ) = ...
                    sum( cm.n( b_todo, cm.b_lambda ) .* reshape( cm.tau( cm.b_lambda ), [ 1, cm.d ] ), 2 );
                
                cm.b_tau_n = cm.b_n;
                
            end
            
        end
        
        function update_p_n ( cm )
            % update cm.p_n
            
            b_todo = cm.b_n & ~cm.b_p_n;
            
            n_todo = sum( b_todo );
            
            if ( n_todo > 0 )
                
                cm.p_n( :, b_todo ) = sum( reshape( cm.n( b_todo, cm.b_lambda ), [ 1, n_todo, cm.d ] ) .* ...
                    reshape( cm.p( :, cm.b_lambda ), [ 3, 1, cm.d ] ), 3 );
                
                cm.b_p_n = cm.b_n;
                
            end
            
        end
        
        function update_En_exp ( cm, b_todo, n_todo, idx_time )
            % update relaxation matrices (in presence of diffusion)
                                    
            if ( cm.diffusion == CoMoTk.Iso_Diff )

                tau_D = cm.tau( idx_time ) .* cm.D;
                
                pn2 = reshape( sum( cm.p_n( :, b_todo ) .* cm.p_n( :, b_todo ), 1 ), [ 1, 1, n_todo ] );
                pns =  reshape( sum( cm.p_n( :, b_todo ) .* cm.s( :, idx_time ), 1 ), [ 1, 1, n_todo ] );

                En_exp_tmp = zeros( 3, 1, n_todo );

                En_exp_tmp( 1, 1, : ) = - ( pn2 + pns + cm.S( idx_time ) );
                En_exp_tmp( 2, 1, : ) = - pn2;
                En_exp_tmp( 3, 1, : ) = - ( pn2 - pns + cm.S( idx_time ) );
                
                % derivatives
                
                for i = 1 : cm.len.dD
                    
                    cm.dEn_exp_dD( :, 1, b_todo, 1, idx_time ) = ...
                        cm.tau( idx_time ) .* En_exp_tmp( :, 1, : );
                        
                end                
                
                if ( cm.len.dtau > 0 )
                    
                    [ ~, i ] = find( cm.dtau == cm.lambda( idx_time ) );
                    
                    if ( ~isempty( i ) )
                        
                        cm.dEn_exp_dtau( :, :, b_todo, 1, i ) = cm.D .* En_exp_tmp( :, 1, : );
                        
                    end
                    
                end

                for i = 1 : cm.len.dp
                    
                    [ ~, i_ ] = find( cm.dp( i ) == cm.lambda );
                    
                    if ( ~isempty( i_ ) )

                        dpn2_dp = 2 .* ...
                            reshape( cm.n( b_todo, i_ ), [ 1, 1, n_todo ] ) .* ...
                            reshape( cm.p_n( :, b_todo )', [ 1, 1, n_todo, 3 ] );
                        
                        dpns_dp = ...
                            reshape( cm.n( b_todo, i_ ), [ 1, 1, n_todo ] ) .* ...
                            reshape( cm.s( :, idx_time ), [ 1, 1, 1, 3 ] );
                        
                        cm.dEn_exp_dp( 1, :, b_todo, :, i, idx_time ) = - tau_D .* ( dpn2_dp + dpns_dp );
                        cm.dEn_exp_dp( 2, :, b_todo, :, i, idx_time ) = - tau_D .* dpn2_dp;
                        cm.dEn_exp_dp( 3, :, b_todo, :, i, idx_time ) = - tau_D .* ( dpn2_dp - dpns_dp );
                        
                    end
                    
                end
                
                if( cm.len.ds > 0 )
                    
                    [ ~, i ] = find( cm.ds == cm.lambda( idx_time ) );
                    
                    if ( ~isempty( i ) )
                        
                        dpns_ds = reshape( cm.p_n( :, b_todo )', [ 1, 1, n_todo, 3 ] );
                        
                        cm.dEn_exp_ds( 1, :, b_todo, :, i ) = - tau_D .* dpns_ds;
                        
                    end
                    
                end
                
                if ( cm.len.dS > 0 )
                    
                    [ ~, i ] = find( cm.dS == cm.lambda( idx_time ) );
                    
                    if ( ~isempty( i ) )

                        cm.dEn_exp_dS( 1, :, 1, 1, i ) = - tau_D;
                        
                    end
                    
                end
                
                % finalize exponent
                
                cm.En_exp( :, :, b_todo, idx_time ) = tau_D .* En_exp_tmp;
                
            elseif ( cm.diffusion == CoMoTk.Tensor_Diff )

                pnDpn = reshape( sum( sum( ...
                    reshape( cm.p_n( :, b_todo ), [ 3, 1, 1, n_todo ] ) .* ...
                    cm.D .* ...
                    reshape( cm.p_n( :, b_todo ), [ 1, 3, 1, n_todo ] ) ...
                    , 1 ), 2 ), [ 1, cm.n_tissues, n_todo ] ); 

                pnDs = reshape( sum( sum( ...
                    reshape( cm.p_n( :, b_todo ), [ 3, 1, 1, n_todo ] ) .* ...
                    cm.D .* ...
                    reshape( cm.s( :, idx_time ), [ 1, 3 ] ) ...
                    , 1 ), 2 ), [ 1, cm.n_tissues, n_todo ] );
                
                DS = reshape( sum( cm.D_vec .* cm.S( :, idx_time ), 1 ), [ 1, cm.n_tissues ] );

                cm.En_exp( 1, :, b_todo, idx_time ) = - ( pnDpn + pnDs + DS );
                cm.En_exp( 2, :, b_todo, idx_time ) = - pnDpn;
                cm.En_exp( 3, :, b_todo, idx_time ) = - ( pnDpn - pnDs + DS );
                
                % derivatives

                if ( cm.len.dD > 0 )
                    
                    dpnDpn_dD = zeros( 1, 1, n_todo, 6 );
                    
                    dpnDpn_dD( 1, 1, :, CoMoTk.Tensor_idx.xx ) = cm.p_n( 1, b_todo ) .* cm.p_n( 1, b_todo );
                    dpnDpn_dD( 1, 1, :, CoMoTk.Tensor_idx.yy ) = cm.p_n( 2, b_todo ) .* cm.p_n( 2, b_todo );
                    dpnDpn_dD( 1, 1, :, CoMoTk.Tensor_idx.zz ) = cm.p_n( 3, b_todo ) .* cm.p_n( 3, b_todo );
                    dpnDpn_dD( 1, 1, :, CoMoTk.Tensor_idx.xy ) = 2 .* cm.p_n( 1, b_todo ) .* cm.p_n( 2, b_todo );
                    dpnDpn_dD( 1, 1, :, CoMoTk.Tensor_idx.xz ) = 2 .* cm.p_n( 1, b_todo ) .* cm.p_n( 3, b_todo );
                    dpnDpn_dD( 1, 1, :, CoMoTk.Tensor_idx.yz ) = 2 .* cm.p_n( 2, b_todo ) .* cm.p_n( 3, b_todo );

                    dpnDs_dD = zeros( 1, 1, n_todo, 6 );
                    
                    dpnDs_dD( 1, 1, :, CoMoTk.Tensor_idx.xx ) = cm.p_n( 1, b_todo ) .* cm.s( 1, idx_time );
                    dpnDs_dD( 1, 1, :, CoMoTk.Tensor_idx.yy ) = cm.p_n( 2, b_todo ) .* cm.s( 2, idx_time );
                    dpnDs_dD( 1, 1, :, CoMoTk.Tensor_idx.zz ) = cm.p_n( 3, b_todo ) .* cm.s( 3, idx_time );
                    dpnDs_dD( 1, 1, :, CoMoTk.Tensor_idx.xy ) = ...
                        cm.p_n( 1, b_todo ) .* cm.s( 2, idx_time ) + cm.p_n( 2, b_todo ) .* cm.s( 1, idx_time );
                    dpnDs_dD( 1, 1, :, CoMoTk.Tensor_idx.xz ) = ...
                        cm.p_n( 1, b_todo ) .* cm.s( 3, idx_time ) + cm.p_n( 3, b_todo ) .* cm.s( 1, idx_time );
                    dpnDs_dD( 1, 1, :, CoMoTk.Tensor_idx.yz ) = ...
                        cm.p_n( 2, b_todo ) .* cm.s( 3, idx_time ) + cm.p_n( 3, b_todo ) .* cm.s( 2, idx_time );

                    dDS_dD = zeros( 1, 1, 1, 6 );
                    
                    dDS_dD( 1, 1, 1, 1 : 3 ) = cm.S( 1 : 3, idx_time );
                    dDS_dD( 1, 1, 1, 4 : 6 ) = 2 .* cm.S( 4 : 6, idx_time );
                    
                    cm.dEn_exp_dD( 1, 1, b_todo, :, idx_time ) = ...
                        - cm.tau( idx_time ) .* ( dpnDpn_dD + dpnDs_dD + dDS_dD );
                    cm.dEn_exp_dD( 2, 1, b_todo, :, idx_time ) = ...
                        - cm.tau( idx_time ) .* dpnDpn_dD;
                    cm.dEn_exp_dD( 3, 1, b_todo, :, idx_time ) = ...
                        - cm.tau( idx_time ) .* ( dpnDpn_dD - dpnDs_dD + dDS_dD );
                    
                end
                
                if ( cm.len.dtau > 0 )
                    
                    [ ~, i ] = find( cm.dtau == cm.lambda( idx_time ) );
                    
                    if ( ~isempty( i ) )
                        
                        cm.dEn_exp_dtau( :, :, b_todo, 1, i ) = cm.En_exp( :, :, b_todo, idx_time );
                        
                    end
                    
                end

                for i = 1 : cm.len.dp
                    
                    [ ~, i_ ] = find( cm.dp( i ) == cm.lambda );
                    
                    if ( ~isempty( i_ ) )
                        
                        dpnDpn_dp = 2 .* reshape( cm.n( b_todo, i_ ), [ 1, 1, n_todo ] ) .* ...
                            reshape( reshape( ...
                            sum( cm.D .* reshape( cm.p_n( :, b_todo ), [ 1, 3, 1, n_todo ] ), 2 ), ...
                            [ 3, cm.n_tissues * n_todo ] )', [ 1, cm.n_tissues, n_todo, 3 ] );
                        
                        dpnDs_dp = reshape( cm.n( b_todo, i_ ), [ 1, 1, n_todo ] ) .* ...
                            reshape( reshape( ...
                            sum( cm.D .* reshape( cm.s( :, idx_time ), [ 1, 3 ] ), 2 ), ...
                            [ 3, cm.n_tissues ] )', [ 1, cm.n_tissues, 1, 3 ] );
                                                
                        cm.dEn_exp_dp( 1, :, b_todo, :, i, idx_time ) = - cm.tau( idx_time ) .* ( dpnDpn_dp + dpnDs_dp );
                        cm.dEn_exp_dp( 2, :, b_todo, :, i, idx_time ) = - cm.tau( idx_time ) .* dpnDpn_dp;
                        cm.dEn_exp_dp( 3, :, b_todo, :, i, idx_time ) = - cm.tau( idx_time ) .* ( dpnDpn_dp - dpnDs_dp );
                        
                    end
                    
                end                

                if( cm.len.ds > 0 )
                    
                    [ ~, i ] = find( cm.ds == cm.lambda( idx_time ) );
                    
                    if ( ~isempty( i ) )
                        
                        dpnDs_ds = ...
                            reshape( reshape( ...
                            sum( cm.D .* reshape( cm.p_n( :, b_todo ), [ 1, 3, 1, n_todo ] ), 2 ), ...
                            [ 3, cm.n_tissues * n_todo ] )', [ 1, cm.n_tissues, n_todo, 3 ] );
                                                
                        cm.dEn_exp_ds( 1, :, b_todo, :, i ) = - cm.tau( idx_time ) .* dpnDs_ds;
                        
                    end
                    
                end

                if ( cm.len.dS > 0 )
                    
                    [ ~, i ] = find( cm.dS == cm.lambda( idx_time ) );
                    
                    if ( ~isempty( i ) )

                        cm.dEn_exp_dS( 1, :, 1, :, i ) = - cm.tau( idx_time ) .* ...
                            reshape( cm.D_vec', [ 1, cm.n_tissues, 1, 6 ] );
                        
                    end
                    
                end
                                
                % finalize exponent
                
                cm.En_exp( :, :, b_todo, idx_time ) = cm.tau( idx_time ) .* cm.En_exp( :, :, b_todo, idx_time );
                                                
            else
                
                error( 'Unexpected value for ''cm.diffusion''.' );
                
            end
                        
        end
        
        function update_relaxation ( cm, idx_time )
            % update relaxation (without and with diffusion)
            
            % transverse and longitudinal relaxation
            
            if ( ~cm.b_E( idx_time ) )
                
                % first call: initialize relaxation
                
                cm.E( 1, :, idx_time ) = exp( - cm.R1 .* cm.tau( idx_time ) );
                cm.E( 2, :, idx_time ) = exp( - cm.R2 .* cm.tau( idx_time ) );
                
                % mark as set
                
                cm.b_E( idx_time ) = true;
                
            end
            
            % diffusion damping
            
            if ( cm.diffusion ~= CoMoTk.No_Diff )
                
                % identify configurations that require updating
                
                b_todo = cm.b_n & ~cm.b_En( :, idx_time );
                
                n_todo = sum( b_todo );
                
                if ( n_todo > 0 )

                    % calculate the exponent for diffusion damping (isotropic and tensor)
                    % and all required derivatives
                    
                    cm.update_En_exp( b_todo, n_todo, idx_time );
                
                    % update diffusion damping factors
                    
                    cm.En( :, :, b_todo, idx_time ) = ...
                        cm.E( [ 2, 1, 2 ], :, idx_time ) .* exp( cm.En_exp( :, :, b_todo, idx_time ) );
                                        
                    % we are done
                    
                    cm.b_En( b_todo, idx_time ) = true;
                    
                end
                                
            end
            
        end
        
        function [ idx_time, b_n_new, b_up_free_time, b_do_free_time ] = pre_time ( cm, lambda )
            % preparations before executing a time interval
            
            if ( cm.options.verbose )
                
                fprintf( 1, 'pre_time: d = %d\n', cm.d );
                fprintf( 1, 'pre_time: n_conf = %d\n', cm.n_conf );
                
            end
            
            % is this a new time interval?
            
            b_time = cm.b_lambda & ( cm.lambda == lambda );
            idx_time = find( b_time, 1, 'first' );
            
            if ( isempty( idx_time ) ) % new time index
                
                if ( cm.alloc_d == cm.d ) % need to allocate more space first
                    
                    cm.alloc_d = cm.alloc_d + cm.options.alloc_d;
                    
                    [ ~, rng_d ] = cm.reallocate_nd;
                    
                    % extend b_time
                    
                    b_time( rng_d ) = false;
                    
                end
                
                % get the first free place
                
                idx_time = find( ~cm.b_lambda, 1, 'first' );
                
                % update logical indices amd set identifier
                
                b_time( idx_time ) = true;
                cm.b_lambda( idx_time ) = true;
                cm.lambda( idx_time ) = lambda;
                
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
            
            %% bookkeeping for lambda == handle_time
            
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
            
            %% bookkeeping for lambda ~= handle_time
            
            % associated logical and actual indices
            
            sz_nd = [ cm.alloc_n, cm.alloc_d ];
            
            b_other_lambda = cm.b_lambda & ~b_time;
            
            b_up_up_o = b_up_free_time & ~cm.b_up_free & b_other_lambda;
            b_up_do_o = b_up_free_time & ~cm.b_do_free & b_other_lambda;
            b_do_up_o = b_do_free_time & ~cm.b_up_free & b_other_lambda;
            b_do_do_o = b_do_free_time & ~cm.b_do_free & b_other_lambda;
            
            b_up_up_f = b_up_free_time & cm.b_up_free & b_other_lambda;
            b_up_do_f = b_up_free_time & cm.b_do_free & b_other_lambda;
            b_do_up_f = b_do_free_time & cm.b_up_free & b_other_lambda;
            b_do_do_f = b_do_free_time & cm.b_do_free & b_other_lambda;
            
            [ r, c ] = find( b_up_up_o );
            
            if ( ~isempty( r ) )
                
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
                
                cm.b_up_free( sub2ind( sz_nd, cm.idx_up_new( r ), c ) ) = true;
                
            end
            
            [ r, c ] = find( b_up_do_f );
            
            if ( ~isempty( r ) )
                                
                cm.b_do_free( sub2ind( sz_nd, cm.idx_up_new( r ), c ) ) = true;
                
            end
            
            [ r, c ] = find( b_do_up_f );
            
            if ( ~isempty( r ) )
                
                cm.b_up_free( sub2ind( sz_nd, cm.idx_do_new( r ), c ) ) = true;
                
            end
            
            [ r, c ] = find( b_do_do_f );
            
            if ( ~isempty( r ) )
                
                cm.b_do_free( sub2ind( sz_nd, cm.idx_do_new( r ), c ) ) = true;
                
            end
            
            % configuration order remains unchanged
            
            cm.n( cm.idx_up_new( b_up_free_time ), b_other_lambda ) = cm.n( b_up_free_time, b_other_lambda );
            cm.n( cm.idx_do_new( b_do_free_time ), b_other_lambda ) = cm.n( b_do_free_time, b_other_lambda );
            
            %% logical indices of new states
            
            cm.b_up_free( cm.idx_up_new( b_up_free_time ), idx_time ) = true;
            cm.b_do_free( cm.idx_up_new( b_up_free_time ), idx_time ) = false;
            
            cm.b_do_free( cm.idx_do_new( b_do_free_time ), idx_time ) = true;
            cm.b_up_free( cm.idx_do_new( b_do_free_time ), idx_time ) = false;
            
            %% now we can update this as well
            
            cm.b_up_free( b_up_free_time, idx_time ) = false;
            cm.b_do_free( b_do_free_time, idx_time ) = false;
            
            %% prepare to update the new occupation
            
            b_n_new = cm.b_n;
            b_n_new( cm.idx_up_new( b_up_free_time ) ) = true;
            b_n_new( cm.idx_do_new( b_do_free_time ) ) = true;
            
        end
        
        function post_time ( cm, b_n_new )
            % several parameter updates after time interval
            
            cm.b_n = b_n_new;
            cm.n_conf = sum( cm.b_n );
            
            if ( cm.options.verbose )
                
                fprintf( 1, 'post_time: d = %d\n', cm.d );
                fprintf( 1, 'post_time: n_conf = %d\n', cm.n_conf );
                
            end
            
        end
                
        function b_is_ok = meltdown ( cm )
            % remove negligible configurations
            
            b_is_ok = true;
            
            % no removal of configurations for cm.epsilon == 0
            
            if ( cm.epsilon == 0 )
                
                return;
                
            end
            
            % calculate size of configuration vectors
            
            abs_m_2 = zeros( cm.alloc_n, 1 );
            abs_m_2( cm.b_n ) = real( sum( reshape( real( cm.m( :, :, cm.b_n ) .* conj( cm.m( :, :, cm.b_n ) ) ), [ 3 * cm.n_tissues, cm.n_conf ] ), 1 ) ).';
            b_small = abs_m_2 < cm.epsilon^2;

            % needed a lot below
            
            sz_nd = [ cm.alloc_n, cm.alloc_d ];
            
            % we do not want to remove these (to be updated later)
            
            b_exclude_n = ~cm.b_n;
            b_exclude_n( cm.null_idx ) = true;
            
            n_candidates_tot = 0;
            n_removed_tot = 0;
            n_cycles = 0;
            
            % removal may require several iterations (from the boundary inwards)
            
            while( true )
                
                %% candidates for removal:
                % (a) occupied
                % (b) not explicitly excluded
                % (c) small
                % (d) have an unoccupied neighbor in every direction
                
                b_remove_n = ...
                    cm.b_n & ...             % (a)
                    ~b_exclude_n & ...       % (b)
                    b_small & ...            % (c)
                    ( sum( cm.b_up_free( :, cm.b_lambda ) | cm.b_do_free( :, cm.b_lambda ), 2 ) == cm.d );  % (d)
                
                n_remove = sum( b_remove_n );
                
                n_candidates_tot = n_candidates_tot + n_remove;
                
                if ( n_remove == 0 )
                    
                    break;  % nothing more to to
                    
                end
                                
                n_removed = 0;
                
                if ( cm.rapid_meltdown ) % looking at the candidates at once (requires more memory)
                
                    % first we have to check that no occupation patterns like
                    %          ---------
                    %          | 1 | 0 |
                    %          ---------
                    %          | 0 | 1 |
                    %          ---------
                    % arise after removal of the configurations

                    idx_n = find( b_remove_n );
                    idx_lambda = find( cm.b_lambda );

                    b_up_occ = ~cm.b_up_free( b_remove_n, cm.b_lambda );
                    b_do_occ = ~cm.b_do_free( b_remove_n, cm.b_lambda );

                    idx_up_occ = cm.idx_up( b_remove_n, cm.b_lambda );
                    idx_do_occ = cm.idx_do( b_remove_n, cm.b_lambda );

                    % up up

                    lambda_pairs = ...
                        reshape( idx_lambda, [ 1, cm.d ] ) > reshape( idx_lambda, [ 1, 1, cm.d ] );
                    
                    [ r_up_up, c_up_up ] = ...
                        find( reshape( ...
                            b_up_occ & reshape( b_up_occ, [ n_remove, 1, cm.d ] ) & lambda_pairs, ...
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
                            reshape( idx_lambda( j_up ), sz_n_ ) ) );
                        
                        b_exclude_n( idx_n( r_up_up( b_free ) ) ) = true;
                        
                        sz_n_ = [ sum( ~b_free ), 1 ];
                        
                        b_exclude_n( cm.idx_up( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_up_occ ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_up_up( ~b_free), ...
                            i_up( ~b_free ) ) ), sz_n_ ), ...
                            reshape( idx_lambda( j_up( ~b_free ) ), sz_n_ ) ) ) ) = true;
                        
                    end
                    
                    % do do
                    
                    [ r_do_do, c_do_do ] = ...
                        find( reshape( ...
                            b_do_occ & reshape( b_do_occ, [ n_remove, 1, cm.d ] ) & lambda_pairs, ...
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
                            reshape( idx_lambda( j_do ), sz_n_ )  ) );
                        
                        b_exclude_n( idx_n( r_do_do( b_free ) ) ) = true;

                        sz_n_ = [ sum( ~b_free ), 1 ];
                                                
                        b_exclude_n( cm.idx_do( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_do_occ ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_do_do( ~b_free), ...
                            i_do( ~b_free ) ) ), sz_n_ ) , ...
                            reshape( idx_lambda( j_do( ~b_free ) ), sz_n_ ) ) ) ) = true;
                        
                    end
                             
                    % up do
                    
                    lambda_pairs = ...
                        reshape( idx_lambda, [ 1, cm.d ] ) ~= reshape( idx_lambda, [ 1, 1, cm.d ] );
                    
                    [ r_up_do, c_up_do ] = ...
                        find( reshape( ...
                            b_up_occ & reshape( b_do_occ, [ n_remove, 1, cm.d ] ) & lambda_pairs, ...
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
                            reshape( idx_lambda( j_do ), sz_n_ )  ) );
                        
                        b_exclude_n( idx_n( r_up_do( b_free ) ) ) = true;
                        
                        sz_n_ = [ sum( ~b_free ), 1 ];
                        
                        b_exclude_n( cm.idx_do( ...
                            sub2ind( sz_nd, ...
                            reshape( idx_up_occ ( ...
                            sub2ind( [ n_remove, cm.d ], ...
                            r_up_do( ~b_free), ...
                            i_up( ~b_free ) ) ), sz_n_ ) , ...
                            reshape( idx_lambda( j_do( ~b_free ) ), sz_n_ ) ) ) ) = true;
                        
                    end
                    
                    % update configurations to be removed
                    
                    b_remove_n( b_exclude_n ) = false;
                    n_removed = sum( b_remove_n );
                    
                    % now we may remove the configurations
                    
                    cm.b_n( b_remove_n ) = false;
                    cm.n_conf = cm.n_conf - n_removed;
                    
                    % disconnect adjacent configurations
                    
                    b_up_occ = b_remove_n & cm.b_lambda & ~cm.b_up_free;
                    b_do_occ = b_remove_n & cm.b_lambda & ~cm.b_do_free;
                    
                    [ r_up_occ, c_up_occ ] = find( b_up_occ );
                    [ r_do_occ, c_do_occ ] = find( b_do_occ );
                    
                    idx_up_occ = cm.idx_up( b_up_occ );
                    idx_do_occ = cm.idx_do( b_do_occ );
                    
                    cm.b_up_free( sub2ind( sz_nd, idx_do_occ, c_do_occ ) ) = true;
                    cm.b_do_free( sub2ind( sz_nd, idx_up_occ, c_up_occ ) ) = true;
                    cm.b_up_free( sub2ind( sz_nd, r_up_occ, c_up_occ ) ) = true;
                    cm.b_do_free( sub2ind( sz_nd, r_do_occ, c_do_occ ) ) = true;
                    
                    cm.b_tau_n( b_remove_n ) = false;
                    cm.b_p_n( b_remove_n ) = false;

                    if ( cm.diffusion ~= CoMoTk.No_Diff )
                        
                        cm.b_En( b_remove_n, cm.b_lambda ) = false;
                        
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
                        
                        b_up_occ = cm.b_lambda & ~cm.b_up_free( r( i ), : );
                        b_do_occ = cm.b_lambda & ~cm.b_do_free( r( i ), : );
                        
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
                        
                        n_removed = n_removed + 1;
                        
                        cm.b_n( r( i ) ) = false;
                        cm.n_conf = cm.n_conf - 1;
                        
                        % disconnect with adjacent configurations
                        
                        cm.b_up_free( sub2ind( sz_nd, idx_do_occ, c_do_occ ) ) = true;
                        cm.b_do_free( sub2ind( sz_nd, idx_up_occ, c_up_occ ) ) = true;
                        cm.b_up_free( r( i ), cm.b_lambda ) = true;
                        cm.b_do_free( r( i ), cm.b_lambda ) = true;
                        
                        cm.b_tau_n( r( i ) ) = false;
                        cm.b_p_n( r( i ) ) = false;
                        
                        if ( cm.diffusion ~= CoMoTk.No_Diff )
                        
                            cm.b_En( r( i ), cm.b_lambda ) = false;
                            
                        end
                                             
                    end
                    
                end
                
                % dimensions without nonzero configurations are in fact
                % nonexistent and can be removed
                
                b_remove_lambda = cm.b_lambda & ~any( cm.n( cm.b_n, : ) ~= 0 );
                
                if ( any( b_remove_lambda ) )
                    
                    if ( cm.options.verbose )
                        
                        fprintf( 1, 'meltdown: Remove %d dimensions.\n', sum( b_remove_lambda ) );
                        
                    end
                    
                    cm.b_lambda( b_remove_lambda ) = false;
                    cm.d = sum( cm.b_lambda );
                    
                    cm.b_up_free( :, b_remove_lambda ) = true;
                    cm.b_do_free( :, b_remove_lambda ) = true;
                    
                    cm.b_tau( b_remove_lambda ) = false;
                    cm.b_p( b_remove_lambda ) = false;
                    
                    cm.b_E( b_remove_lambda ) = false;
                    
                    if ( cm.diffusion ~= CoMoTk.No_Diff )
                        
                        cm.b_En( :, b_remove_lambda ) = false;
                        
                    end
                    
                end                
                
                n_cycles = n_cycles + 1;
                n_removed_tot = n_removed_tot + n_removed;
                
            end
                        
            if ( cm.options.verbose )
                    
                fprintf( 1, 'meltdown: removed %s of %s candidates in %s cycles.\n', num2str( n_removed_tot ), num2str( n_candidates_tot ), num2str( n_cycles ) );
                    
            end
                
        end
        
        function [ rng_n, rng_d ] = reallocate_nd ( cm )
            % reallocate new space, if necessary
            
            [ alloc_n_old, alloc_d_old ] = size( cm.n );
            
            if ( cm.options.debug )
                
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
                
                if ( cm.diffusion ~= CoMoTk.No_Diff )
                    
                    cm.b_En( rng_n, : ) = false;
                    cm.En( :, :, rng_n, : ) = 0;
                    cm.En_exp( :, :, rng_n, : ) = 0;
                    
                end
                
                cm.b_tau_n( rng_n ) = false;
                cm.tau_n( rng_n ) = 0;
                
                cm.b_p_n( rng_n ) = false;
                cm.p_n( :, rng_n ) = 0;
                
                if ( cm.len.dR1 > 0 )
                    
                    cm.dm.dR1( :, :, rng_n, 1, : ) = 0;
                    
                end
                
                if ( cm.len.dR2 > 0 )
                    
                    cm.dm.dR2( :, :, rng_n, 1, : ) = 0;
                    
                end
                
                if ( cm.len.dD > 0 )
                    
                    cm.dm.dD( :, :, rng_n, :, : ) = 0;
                    cm.dEn_exp_dD( :, :, rng_n, :, : ) = 0;
                    
                end
                
                if ( cm.len.dB1 > 0 )
                    
                    cm.dm.dB1( :, :, rng_n ) = 0;
                    
                end
                
                if ( cm.len.dFlipAngle > 0 )
                    
                    cm.dm.dFlipAngle( :, :, rng_n, 1, : ) = 0;
                    
                end
                
                if ( cm.len.dPhase > 0 )
                    
                    cm.dm.dPhase( :, :, rng_n, 1, : ) = 0;
                    
                end
                
                if ( cm.len.dtau > 0 )
                    
                    cm.dm.dtau( :, :, rng_n, 1, : ) = 0;
                    cm.dEn_exp_dtau( :, :, rng_n, 1, : ) = 0;
                    
                end
                
                if ( cm.len.dp > 0 )
                    
                    cm.dm.dp( :, :, rng_n, :, : ) = 0;
                    cm.dEn_exp_dp( :, :, rng_n, :, :, : ) = 0;
                    
                end
                
                if ( cm.len.ds > 0 )
                    
                    cm.dm.ds( :, :, rng_n, :, : ) = 0;
                    cm.dEn_exp_ds( :, :, rng_n, :, : ) = 0;
                    
                end
                
                if ( cm.len.dS > 0 )
                    
                    cm.dm.dS( :, :, rng_n, :, : ) = 0;
                    cm.dEn_exp_dS( :, :, rng_n, :, : ) = 0;
                    
                end
                
            end
            
            if ( ~isempty( rng_d ) )
                
                cm.b_lambda( rng_d ) = false;
                cm.lambda( rng_d ) = 0;
                
                cm.n( :, rng_d ) = 0;
                
                cm.b_up_free( :, rng_d ) = true;
                cm.b_do_free( :, rng_d ) = true;
                
                cm.idx_up( :, rng_d ) = 0;
                cm.idx_do( :, rng_d ) = 0;
                
                cm.b_E( rng_d ) = false;
                cm.E( :, :, rng_d ) = 0;
                
                if ( cm.diffusion ~= CoMoTk.No_Diff )
                    
                    cm.b_En( :, rng_d ) = false;
                    cm.En( :, :, :, rng_d ) = 0;
                    cm.En_exp( :, :, :, rng_d ) = 0;
                    
                    cm.s( :, rng_d ) = 0;
                    cm.S( :, rng_d ) = 0;
                    
                end

                cm.b_tau( rng_d ) = false;
                cm.tau( rng_d ) = 0;
                
                cm.b_p( rng_d ) = false;
                cm.p( :, rng_d ) = 0;
                                
                if ( cm.len.dD > 0 )
                    
                    cm.dEn_exp_dD( :, :, :, :, rng_d ) = 0;
                    
                end
                
                if ( cm.len.dp > 0 )
                    
                    cm.dEn_exp_dp( :, :, :, :, :, rng_d ) = 0;
                    
                end
                
            end
            
            if ( cm.options.debug )
                
                fprintf( 1, 'reallocate_nd: (after) alloc_n = %d\n', cm.alloc_n );
                fprintf( 1, 'reallocate_nd: (after) alloc_d = %d\n', cm.alloc_d );
                
            end
            
        end
        
    end
    
end