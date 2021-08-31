classdef CoMo < matlab.mixin.Copyable
    % CoMo  Configuration Model
    %
    % Implementation of the configuration model (CM), which has
    % been extensively described in a separate manuscript.
    % (link available on: https://github.com/cganter/CoMo)
    % 
    % Also check out the specific documentation and example scripts.
    %
    % Requires Matlab R2018b (or later).
    %
    % (c) Carl Ganter 2020
    
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
        
        % default preallocation
        
        alloc = 1000;
        
        % default accuracy (0 == perfect)
        
        epsilon = 0;
        
        % some informative output (occupied space, resizing, reallocations)
        
        verbose = false;
        
    end
    
    properties ( Dependent )
        
        % ALL dependent quantities (except for the optional s)
        % have to be set explicitly prior to simulation. Otherwise an error is returned.
        
        %% tissue properties
        % relaxation rates
        
        R1;                                                                % [ 1 / ms ]
        R2;                                                                % [ 1 / ms ]
        
        % diffusion sensitivity (scalar or tensor)
        
        D;                                                                 % [ um^2 / ms ]  i.e. D( Water ) is approximately equal to 3
        
        %% multiple tissues
        % (relative) proton density
        
        pd;
        
        % delta omega (due to chemical shift, e.g. to specify water/fat models)
        
        dom;                                                               % [ rad / ms ]
        
        % transition rates (magnetization transfer and exchange)
        
        k;                                                                 % [ 1 / ms ]
        
        % as unit for bulk motion we assume [ mm / s ]
        % (used in method "time")
        
        %% measurement conditions
        % relative B1
        
        B1;
        
        %% discretization of (p,tau)-space
        
        d_p;
        d_tau;
                
    end
    
    properties ( SetAccess = private )
        
        b_prep = false;                                                    % is everything prepared?
        
        %% information about the actual state
        
        % The configuration vectors are stored in a conventional 2d
        % matrix. Each of the two dimensions merges several subdimensions in the following order:
        %
        % dimension 1: 3d components (spatial and gradient moments), tissue
        % dimension 2: tau, p (together 4 components)
        
        % stores information, which part of the allocated space is actually
        % occupied by configurations, type == bool, 
        % size( cm.occ ) = [ 1, cm.alloc ]
        
        occ = [];
        
        % associated information

        n_occ = [];              % number of stored configurations
        occ_idx = [];            % location in allocated space == find( cm.occ )
        null_idx = [];           % location, related to the origin of configuration space (0,0)
                                 % (always occupied, regardless of the configuration value)        
        
        % each occupied configuration is specified by the pair ( tau, p )
        % the array cm.occ_sub of size( cm.occ_sub ) == [ 4, cm.alloc ]
        % stores integers, such that, for example, the actual value of tau
        % is given by
        % cm.tau = cm.d_tau * cm.occ_sub( :, 1 )
        
        occ_sub = [];
        
        % at any time, these infos about the occupied configurations are 
        % guaranteed to be consistent
                
        m = [];                                                            % configuration vector
                                                                           % m = 6d-complex array
                                                                           % dimensions:
                                                                           % 1: magnetization vector components (== 3)
                                                                           % 2: number of tissues (== cm.n_tissues)
                                                                           % 3-5: accumulated moments p(j) (== 2 * cm.n_p(j) + 1)
                                                                           % 6: accumulated transverse time tau (== 2 * cm.n_tau + 1)
                                                                           
        m_dim = [];                                                        % size( cm.m )
        
        n_tissues = 0;		                                               % number of tissues

        p = [];                                                            % accumulated gradient moments in configuration space
                                                                           % size( cm.p ) == [ 3, cm.alloc ]
        
        tau = [];                                                          % accumulated net transverse durations
                                                                           % size( cm.tau ) == [ 1, cm.alloc ]
    
        %% mandatory tissue properties
        % general properties
        
        diffusion = CoMo.No_Diff;
        
        %% associated with dependent variables
        
        R1_priv = [];
        R2_priv = [];
        D_priv = [];
           
        pd_priv = 1;
        dom_priv = 0;
        k_priv = [];
        B1_priv = 1;

        d_p_priv = zeros( 1, 3 );                                          % 3-d real, associated small increments (absolute values)
        n_p_priv = zeros( 1, 3 );                                          % 3-d integer, number of positive gradient moments in each direction
        
        d_tau_priv = 0;                                                    % real, associated small increment (absolute value)
        n_tau_priv = 0;                                                    % integer, number of positive accumulated times

        %% auxiliary variables
        
        % tensor diffusion
        
        D_vec = [];
        
        % some information is logged here
        
        log = struct( ...
            't', [], ...                                                   % evolved time
            'n_occ', [], ...                                               % number of occupied configurations
            'n_del', [] ...                                                % deleted configurations (for finite epsilon)
        );
        
    end
    
    methods
        
        %% constructor
        
        function cm = CoMo ( )
            % constructor (no arguments, no action)
            
        end
        
        %% set / get methods (mandatory variables)
        
        function set.R1 ( cm, R1 )
            % set relaxation rate R1 = 1 / T1
            %
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMo.prep.' );
                
            end
            
            if ( min( R1 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R1 );
                cm.pd = ones( 1, cm.n_tissues );
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
                
                error( 'Tissue parameters can only be set prior to calling CoMo.prep.' );
                
            end
            
            if ( min( R2 ) < 0 )
                
                error( 'Relaxation rate must not be negative.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( R2 );
                cm.pd = ones( 1, cm.n_tissues );
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
                
                error( 'Tissue parameters can only be set prior to calling CoMo.prep.' );
                
            end
            
            sz_D = size( D );
            
            if ( sz_D( 1 ) ~= 3 || sz_D( 2 ) ~= 3 )
                
                cm.diffusion = CoMo.Iso_Diff;
                
                n_tis = sz_D( 2 );
                
            else
                
                cm.diffusion = CoMo.Tensor_Diff;
                
                if ( length( sz_D ) == 3 )
                    
                    n_tis = sz_D( 3 );
                    
                else
                    
                    n_tis = 1;
                    
                end
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = n_tis;
                cm.pd = ones( 1, cm.n_tissues );
                cm.dom = zeros( cm.n_tissues, 1 );
                
            elseif ( n_tis ~= cm.n_tissues )
                
                error( 'Supplied D does not match number of tissues.' );
                
            end
            
            if ( cm.diffusion == CoMo.Tensor_Diff )
                
                cm.D_vec = zeros( 6, cm.n_tissues );
                
                cm.D_vec( CoMo.Tensor_idx.xx, : ) = D( 1, 1, : );
                cm.D_vec( CoMo.Tensor_idx.yy, : ) = D( 2, 2, : );
                cm.D_vec( CoMo.Tensor_idx.zz, : ) = D( 3, 3, : );
                cm.D_vec( CoMo.Tensor_idx.xy, : ) = 2 .* D( 1, 2, : );
                cm.D_vec( CoMo.Tensor_idx.xz, : ) = 2 .* D( 1, 3, : );
                cm.D_vec( CoMo.Tensor_idx.yz, : ) = 2 .* D( 2, 3, : );
                
            end
            
            cm.D_priv = D;
            
            if ( max( abs( cm.D( : ) ) ) == 0 )
                
                cm.diffusion = CoMo.No_Diff;
                
            end
            
        end
        
        function res = get.D ( cm )
            % return diffusion constant or tensor
            
            res = cm.D_priv;
            
        end
        
        function set.d_p ( cm, d_p )
            % set positive gradient moment increment in each direction
            %
            % IN
            %
            % non-negative vector of length 3 
            
            %             if ( cm.b_prep )
            %
            %                 error( 'Tissue parameters can only be set prior to calling CoMo.prep.' );
            %
            %             end
            
            if ( length( d_p ) ~= 3 || min( d_p ) < 0 )
                
                error( 'Invalid setting of d_p.' );
                
            end
                        
            cm.d_p_priv = d_p( : )';
            
        end
        
        function res = get.d_p ( cm )
            
            res = cm.d_p_priv;
            
        end        
        
        function set.d_tau ( cm, d_tau )
            % set positive gradient moment increment in each direction
            %
            % IN
            %
            % positive number 
            
            %             if ( cm.b_prep )
            %
            %                 error( 'Tissue parameters can only be set prior to calling CoMo.prep.' );
            %
            %             end
            
            if ( d_tau <= 0 )
                
                error( 'Invalid setting of d_tau.' );
                
            end
                        
            cm.d_tau_priv = d_tau;
            
        end
        
        function res = get.d_tau ( cm )
            
            res = cm.d_tau_priv;
            
        end        
           
        %% set / get methods (optional variables)
        
        % relative B1+
        
        function set.B1 ( cm, B1 )
            % set relative B1+
            %
            % set to 1, if not explicitly specified            
            
            if ( cm.b_prep )
                
                error( 'B1 can only be set prior to calling CoMo.prep.' );
                
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
        
        function set.pd ( cm, pd )
            % set the (relative) proton density
            %
            % default = 1 (only, if cm.n_tissues == 1)
            % mandatory parameter, if cm.n_tissues > 1
            % length must match the number of subtissues (e.g. for fat models)
            
            if ( cm.b_prep )
                
                error( 'Tissue parameters can only be set prior to calling CoMo.prep.' );
                
            end
            
            if ( cm.n_tissues == 0 )
                
                cm.n_tissues = length( pd );
                
            elseif ( length( pd ) ~= cm.n_tissues )
                
                error( 'length( pd ) must be equal to number of tissues.' );
                
            end
            
            cm.pd_priv = reshape( pd, [ 1, cm.n_tissues ] );
            
        end
        
        function res = get.pd ( cm )
            % return the (relative) proton density
            
            res = cm.pd_priv;
            
        end
        
        % frequency offset due to chemical shift
        
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
                
        % transition rate matrix (ME/MT)
        
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
            
            % calculate the rotation matrix
            
            c_al = cos( ActualFlipAngle );
            s_al = sin( ActualFlipAngle );
            
            ei_ph = cos( param.Phase ) + 1i * sin( param.Phase );
            
            RotMat = zeros( 3 );
            
            RotMat( 1, 1 ) = 0.5 * ( 1 + c_al );
            RotMat( 1, 2 ) = - 1i * s_al * ei_ph * CoMo.sqrt_0p5;
            RotMat( 1, 3 ) = 0.5 * ( 1 - c_al ) * ei_ph^2;
            RotMat( 2, 1 ) = - conj( RotMat( 1, 2 ) );
            RotMat( 2, 2 ) = c_al;
            RotMat( 2, 3 ) = - RotMat( 1, 2 );
            RotMat( 3, 1 ) = conj( RotMat( 1, 3 ) );
            RotMat( 3, 2 ) = - RotMat( 2, 1 );
            RotMat( 3, 3 ) = RotMat( 1, 1 );

            %% Magnetization exchange/transfer
            
            b_MT = false( 1, cm.n_tissues );
            n_MT = 0;
            
            if ( cm.n_tissues > 1 )
                 
                % check for magnetization transfer (MT)
                
                if ( isfield( param, 'MT_sat' ) )
                    
                    if ( length( param.MT_sat ) ~= cm.n_tissues )
                        
                        error( 'wrong size of saturation due to magnetization transfer' );
                        
                    end
                    
                    % a value > 0 is treated as a bound pool
                    
                    b_MT( : ) = param.MT_sat > 0;
                    n_MT = sum( b_MT );

                end
                
                cRM = cell( cm.n_tissues - n_MT, 1 );
                cRM( : ) = { sparse( RotMat ) };
                
                RotMat = blkdiag( cRM{ : } );
                
            end
                
            % check for magnetization transfer (MT)
            
            if ( n_MT > 0 )
                
                % to access the correct indices
                
                b_ind_MT = repmat( b_MT, [ 3, 1 ] );
                
                if ( cm.n_tissues > n_MT ) % should be always the case..
                
                    cm.m( ~b_ind_MT( : ), cm.occ ) = ...
                        RotMat * cm.m( ~b_ind_MT( : ), cm.occ );
                    
                end
 
                % MT saturation affects only the longitudinal component

                b_MT_long = b_ind_MT & repmat( [ false; true; false ], [ 1, cm.n_tissues ] );
                
                cm.m( b_MT_long( : ), cm.occ ) = ...
                    spdiags( param.MT_sat( b_MT ), 0, n_MT, n_MT ) * ...
                    cm.m( b_MT_long( : ), cm.occ );

                %                 cm.m( b_ind_MT( 2, : ), cm.occ ) = ...
                %                     spdiags( param.MT_sat( b_MT ), 0, n_MT, n_MT ) * ...
                %                     cm.m( b_ind_MT( 2, : ), cm.occ );
                
                % since transverse magnetization is never generated, we
                % don't need to set it to zero explicitly

            else
                
                %% update the configuration cm.m
                
                cm.m( :, cm.occ ) = RotMat * cm.m( :, cm.occ );
                
            end
                               
        end
        
        function time ( cm, param )
            % Applies time interval
            %
            % IN
            %
            % param : structure with the info about the RF pulse,
            % specifically:
            %
            % param.tau    : duration of interval [ ms ]
            % (mandatory parameter)
            % 
            % param.p      : zero-order gradient moment [ rad / um ]
            % (set to zero, if absent)
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
            % (optional, only relevant in case of bulk motion)za

            %% prepare everything, if not done so yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            %% check for duration and zero-order gradient moment
            
            % duration (mandatory)
            
            if ( isfield( param, 'tau' ) )
                
                tau_ = param.tau;
                
            else
                
                error( 'Time interval requires duration.' );
                
            end
            
            % zero order gradient moment
            
            if ( isfield( param, 'p' ) )
                
                p_ = param.p( : );
                
            else
                
                p_ = zeros( 3, 1 );
                
            end
            
            %% prepare relaxation
            
            E = zeros( 3, cm.n_tissues );
            
            E( 1, : ) = exp( - cm.R2 .* tau_ );
            E( 2, : ) = exp( - cm.R1 .* tau_ );
            E( 3, : ) = E( 1, : );
            
            %% prepare diffusion damping
            
            % first, set gradient shape
            
            if ( cm.diffusion ~= CoMo.No_Diff )
                
                if ( isfield( param, 's' ) )
                                        
                    s_ = param.s( : );
                    
                else
                    
                    s_ = p_;
                    
                end
                
                if ( isfield( param, 'S' ) )
                    
                    if ( cm.diffusion == CoMo.Iso_Diff )
                        
                        S_vec = param.S;
                        
                    elseif ( cm.diffusion == CoMo.Tensor_Diff )
                        
                        S_vec = zeros( 6, 1 );
                        
                        S_vec( CoMo.Tensor_idx.xx ) = param.S( 1, 1 );
                        S_vec( CoMo.Tensor_idx.yy ) = param.S( 2, 2 );
                        S_vec( CoMo.Tensor_idx.zz ) = param.S( 3, 3 );
                        S_vec( CoMo.Tensor_idx.xy ) = param.S( 1, 2 );
                        S_vec( CoMo.Tensor_idx.xz ) = param.S( 1, 3 );
                        S_vec( CoMo.Tensor_idx.yz ) = param.S( 2, 3 );
                        
                    else
                        
                        error( 'This should not happen.' );
                        
                    end
                                        
                else  % assume constant gradient during time interval
                    
                    if ( cm.diffusion == CoMo.Iso_Diff )
                        
                        S_vec = ( p_' * p_ ) / 3;
                        
                    elseif ( cm.diffusion == CoMo.Tensor_Diff )
                        
                        S_vec = zeros( 6, 1 );
                        
                        S_tmp = ( p_ * p_' ) ./ 3;
                        
                        S_vec( CoMo.Tensor_idx.xx ) = S_tmp( 1, 1 );
                        S_vec( CoMo.Tensor_idx.yy ) = S_tmp( 2, 2 );
                        S_vec( CoMo.Tensor_idx.zz ) = S_tmp( 3, 3 );
                        S_vec( CoMo.Tensor_idx.xy ) = S_tmp( 1, 2 );
                        S_vec( CoMo.Tensor_idx.xz ) = S_tmp( 1, 3 );
                        S_vec( CoMo.Tensor_idx.yz ) = S_tmp( 2, 3 );
                        
                    else
                        
                        error( 'This should not happen.' );
                        
                    end
                    
                end
                
                S_ = S_vec( : );

            end
            
            % now calculate the damping factors

            Diff_Damp = [];    % if not empty: size( Diff_Damp ) == cm.m_dim
            
            if ( cm.diffusion ~= CoMo.No_Diff )
                
                if ( cm.diffusion == CoMo.Iso_Diff )
                    
                    pn2 = sum( cm.p( :, cm.occ ) .* cm.p( :, cm.occ ) );   % size( pn2 ) == [ 1, cm.n_occ ]
                    pns = sparse( s_' ) * cm.p( :, cm.occ );               % size( pns ) == [ 1, cm.n_occ ]
                    
                    tau_D = tau_ .* cm.D( : );                             % size( tau_D ) == [ cm.n_tissues, 1 ]

                    Diff_Damp = zeros( cm.m_dim, cm.n_occ );
                    
                    Diff_Damp( [ 1 : 3 : end, 2 : 3 : end, 3 : 3 : end ], : ) = ...
                        exp( [ ...
                        - tau_D .* ( pn2 + pns + S_ ); ...
                        - tau_D .* pn2; ...
                        - tau_D .* ( pn2 - pns + S_ ); ...
                        ] );
                    
                elseif ( cm.diffusion == CoMo.Tensor_Diff )
                    
                    tau_D = tau_ .* cm.D_vec';                             % size == [ cm.n_tissues, 6 ];
                    
                    pp = [ ...                                             % size == [ 6, cm.n_occ ]
                        cm.p( 1, cm.occ ) .* cm.p( 1, cm.occ ); ...        % zz
                        cm.p( 2, cm.occ ) .* cm.p( 2, cm.occ ); ...        % yy
                        cm.p( 3, cm.occ ) .* cm.p( 3, cm.occ ); ...        % zz
                        cm.p( 1, cm.occ ) .* cm.p( 2, cm.occ ); ...        % xy
                        cm.p( 1, cm.occ ) .* cm.p( 3, cm.occ ); ...        % xz
                        cm.p( 2, cm.occ ) .* cm.p( 3, cm.occ ) ];          % yz
                    
                    tau_pDp = tau_D * pp;                                  % size == [ cm.n_tissues, cm.n_occ ]
                    
                    tau_D = tau_ .* ...
                        reshape( s_' * reshape( cm.D, 3, [] ), 3, [] )';   % size == [ cm.n_tissues, 3 ];
                    
                    tau_pDs = tau_D * cm.p( :, cm.occ );                   % size == [ cm.n_tissues, cm.n_occ ] 
                    
                    tau_DS_ = tau_ .* ( cm.D_vec' * S_ );                  % size == [ cm.n_tissues, 1 ]
                    
                    Diff_Damp = zeros( cm.m_dim, cm.n_occ );
                    
                    Diff_Damp( [ 1 : 3 : end, 2 : 3 : end, 3 : 3 : end ] , : ) = ...
                        exp( [ ...
                        - tau_pDp - tau_pDs - tau_DS_; ...
                        - tau_pDp; ...
                        - tau_pDp + tau_pDs - tau_DS_ ...
                        ] );
                                                            
                else
                    
                    error( 'Unexpected value for ''cm.diffusion''.' );
                    
                end
                    
            end
            
            %% initialize chemical shift
            
            cs_tau = ones( 3, cm.n_tissues );
            
            cs_tau( 1, : ) = exp( - 1i .* cm.dom .* tau_ );
            cs_tau( 3, : ) = conj( cs_tau( 1, : ) );
            
            %% execute: relaxation, chemical shift, diffusion damping, chemical exchange
            
            if ( isempty( cm.k ) )              % no magnetization exchange
                
                %% update the configuration cm.m
                
                % chemical shift and relaxation                
                
                cm.m( :, cm.occ ) = ( cs_tau( : ) .* E( : ) ) .* cm.m( :, cm.occ );
                
                % diffusion damping, if relevant
                
                if ( ~isempty( Diff_Damp ) )

                    cm.m( :, cm.occ ) = Diff_Damp .* cm.m( :, cm.occ );
                    
                end
                
            else                              % with magnetization exchange

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

                exp_Y_tau_1 = expm( Y_1 .* tau_ );
                exp_Y_tau_0 = expm( Y_0 .* tau_ );

                y = cm.R1( : ) .* cm.pd( : );
                
                % relaxation
                  
                cm.m( 1 : 3 : end, cm.occ ) = ...
                    exp_Y_tau_1 * cm.m( 1 : 3 : end, cm.occ );
                
                cm.m( 2 : 3 : end, cm.occ ) = ...
                    exp_Y_tau_0 * cm.m( 2 : 3 : end, cm.occ );
                
                cm.m( 3 : 3 : end, cm.occ ) = ...
                    conj( exp_Y_tau_1 ) * cm.m( 3 : 3 : end, cm.occ );
                
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
                    
                    dth = - ( 0.5 * tau_ ) * ( p_' * param.v( : ) );
                    
                end
                
                dth = dth - p_' * param.x( : );
                
                ei_dth = exp( 1i * dth );
                
                cm.m( 1 : 3 : end, cm.occ ) = ...
                    ei_dth .* cm.m( 1 : 3 : end, cm.occ );
                
                cm.m( 3 : 3 : end, cm.occ ) = ...
                    conj( ei_dth ) .* cm.m( 3 : 3 : end, cm.occ );
                
            end
            
            %% execute shift in configuration space
            
            % we do this for each direction
            
            tau_shift = round( tau_ / cm.d_tau );
            
            p_shift = zeros( 1, 3 );
            
            for i = 1 : 3
                
                if ( p_( i ) ~= 0 )
                    
                    if ( cm.d_p( i ) == 0 )
                        
                        error( 'Nonzero moment in unexpected direction.' );
                        
                    end
                    
                    p_shift( i ) = round( p_( i ) / cm.d_p( i ) );
                    
                end

            end
            
            % reindex, if necessary            
            
            shift_ = [ tau_shift, p_shift ];

            % calculate shifted occupied indices
            % (in different directions for the two transverse components)
            
            occ_sub_p = cm.occ_sub( cm.occ, : ) + shift_;
            occ_sub_m = cm.occ_sub( cm.occ, : ) - shift_;
            
            % concatenate the occupied (z) and shifted (p,m) locations in
            % configuration space 
            % (i.e. with integer based indexing of tau and p)
            
            os = [ cm.occ_sub( cm.occ, : ); occ_sub_m; occ_sub_p ];

            rng_z = 1 : cm.n_occ;
            rng_m = cm.n_occ + 1 : 2 * cm.n_occ;
            rng_p = 2 * cm.n_occ + 1 : 3 * cm.n_occ;
            
            % sort the rows hierarchically:
            % n-th criterion == value of n-th column            
            
            [ s_os, i_os ] = sortrows( os );
            
            % now we look for the corresponding allocated locations
            % when we are finished, each row in os_idx should map to the
            % allocated location, such that:
            % - the allocated locations, corresponding to cm.occ are not changed
            % - locations from occ_sub_p and occ_sub_m, which were already
            % occupied are assigned the correct allocated locations
            % - the rest is assigned a new location index
            
            os_idx = zeros( 3 * cm.n_occ, 1 );    % there should be no zeros left once we are finished
            os_idx( rng_z ) = cm.occ_idx;         % already occupied locations remain there (for the moment)

            % we now search for multiple entries in configuration space
            % (essentially those shifted locations, which were already occupied)
                        
            % search for rows with equal predecessor
            
            b_eq = [ false; ~any( s_os( 2 : end, : ) - s_os( 1 : end - 1, : ), 2 ) ];
            
            % check that no such row is associated with (z)
            % (just to be sure - can be commented out, once we are)
            
            if ( any( i_os( b_eq ) <= cm.n_occ ) )
                
                error( 'unexpected sortrows behavior - type I' );
                
            end
                
            % these locations are equal to exactly two predecessors
            
            b_eq_2 = [ false; false; b_eq( 3 : end ) & b_eq( 2 : end - 1 ) ];
            
            % while for those only the direct predecessor is equal
            
            b_eq_1 = b_eq & ~b_eq_2;
            
            % the associated locations
            
            i_eq_1 = find( b_eq_1 );
            i_eq_2 = find( b_eq_2 );
            
            % the reference of locations for b_eq_2 should be occupied, if
            % the sorting is not screwed up
            % we confirm that this is not the case
            % (just to be sure - can be commented out, once we are)
            
            %             if ( any( i_os( i_eq_2 - 2 ) > cm.n_occ ) )
            %
            %                 error( 'unexpected sortrows behavior - type II' );
            %
            %             end
            
            % assign the correct occupied index
            
            os_idx( i_os( i_eq_2 ) ) = os_idx( i_os( i_eq_2 - 2 ) );
                
            % the reference for locations in b_eq_1 will also usually be
            % occupied.
            % only in improbable cases of interior holes, generated by discarding of
            % configurations, this might not be the case.
            % to be on the safe side, we therefore pretend this can happen..
            
            b_hole = i_os( i_eq_1 - 1 ) > cm.n_occ;
            
            % if the reference is occupied, we can assign its index
            
            os_idx( i_os( i_eq_1( ~b_hole ) ) ) = os_idx( i_os( i_eq_1( ~b_hole ) - 1 ) );
            
            % determine the required number of new occupied indices.
            % should be equal to the number of remaining zero indices minus the
            % number of holes.
            
            n_new = sum( os_idx == 0 ) - sum( b_hole );
            
            % reallocate, if necessary
            
            if ( n_new + cm.n_occ > cm.alloc )
                
                new_alloc = 2 * cm.alloc;
                
                while ( n_new + cm.n_occ > new_alloc )
                    
                    new_alloc = 2 * new_alloc;
                    
                end
                
                if ( cm.verbose )
                    
                    fprintf( 1, 'Reallocate: %d \t->\t %d\n', cm.alloc, new_alloc );
                    
                end
                
                % resize, but keep indices
                
                occ_tmp = cm.occ;
                cm.occ = false( 1, new_alloc );
                cm.occ( 1 : cm.alloc ) = occ_tmp;
                
                occ_sub_tmp = cm.occ_sub;
                cm.occ_sub = zeros( new_alloc, 4 );
                cm.occ_sub( 1 : cm.alloc, : ) = occ_sub_tmp;
                
                m_tmp = cm.m;
                cm.m = zeros( cm.m_dim, new_alloc );
                cm.m( :, 1 : cm.alloc ) = m_tmp;
                
                tau_tmp = cm.tau;
                cm.tau = zeros( 1, new_alloc );
                cm.tau( :, 1 : cm.alloc ) = tau_tmp;
                
                p_tmp = cm.p;
                cm.p = zeros( 3, new_alloc );
                cm.p( :, 1 : cm.alloc ) = p_tmp;
                
                cm.alloc = new_alloc;
                
            end
            
            % where to store the new configurations?
            
            idx_new = find( ~cm.occ, n_new );

            % now we assign the remaining indices

            % first we create a mask of those locations in b_eq_1, corresponding to
            % holes (see above), since the assignment must be postponed
            
            b_eq_1_hole = b_eq_1;
            b_eq_1_hole( i_eq_1( ~b_hole ) ) = false;
            
            t_ = true( 3 * cm.n_occ, 1 );
            t_( i_os( b_eq_1_hole ) ) = false;
            
            os_idx( os_idx == 0 & t_ ) = idx_new;
            
            % finally take care of assigning indices from new_idx to the holes
            
            os_idx( i_os( i_eq_1( b_hole ) ) ) = os_idx( i_os( i_eq_1( b_hole ) - 1 ) );
            
            % the assignment of indices to os_idx is now finished.
            % it comprises 3 * cm.alloc entries, of which only 
            % cm.alloc + n_new are unique
            % later, these new indices form the entries of 
            % cm.occ_idx = find( cm.occ )
            % since sorting is more time efficient, if the input is already
            % (partly) sorted, we want to maintain an ordered state of
            % cm.occ_sub
            
            map_os_idx = zeros( cm.occ_idx( end ) + n_new, 1 );
            map_os_idx( unique( os_idx( i_os ), 'stable' ) ) = unique( os_idx );
            new_os_idx = map_os_idx( os_idx );
            
            % now perform the shift for the magnetization density
            
            m_ = cm.m( : , cm.occ );
            cm.m( :, cm.occ ) = 0;

            cm.m( 1 : 3 : end, new_os_idx( rng_p ) ) = m_( 1 : 3 : end, : );
            cm.m( 2 : 3 : end, new_os_idx( rng_z ) ) = m_( 2 : 3 : end, : );
            cm.m( 3 : 3 : end, new_os_idx( rng_m ) ) = m_( 3 : 3 : end, : );
            
            % update occupation and indices
            
            cm.n_occ = cm.n_occ + n_new;
            cm.occ( idx_new ) = true;
            cm.occ_idx = find( cm.occ );
            cm.occ_sub( cm.occ, : ) = s_os( ~b_eq, : );
            cm.null_idx = cm.occ_idx( find( ~any( cm.occ_sub( cm.occ, : ), 2 ) ) );
            
            % update tau and p
            
            cm.tau( cm.occ_idx ) = cm.d_tau .* ( cm.occ_sub( cm.occ_idx, 1 ) - cm.occ_sub( cm.null_idx, 1 ) );
            
            for i = 1 : 3
            
                cm.p( i, cm.occ_idx ) = cm.d_p( i ) .* ( cm.occ_sub( cm.occ_idx, i + 1 ) - cm.occ_sub( cm.null_idx, i + 1 ) );
                
            end

            % update rest
            
            %% repolarization
            
            if ( isempty( cm.k ) )              % no magnetization exchange
                                
                cm.m( 2 : 3 : end, cm.null_idx ) = ...
                    cm.m( 2 : 3 : end, cm.null_idx ) + ...
                    ( cm.pd .* ( 1 - E( 2, : ) ) )';
                
            else                              % with magnetization exchange

                % (here, the relative proton densities come into play)
                
                cm.m( 2 : 3 : end, cm.null_idx ) = ...
                    cm.m( 2 : 3 : end, cm.null_idx ) - ...
                    reshape( Y_0 \ ( ( eye( cm.n_tissues ) - exp_Y_tau_0 ) * y ), [ cm.n_tissues, 1 ] );
                
            end            
            
            % now we remove small configurations
            % for epsilon == 0, this just frees locations, which are no
            % longer occupied after shifting
            
            occ_idx_ = cm.occ_idx;
            occ_idx_( occ_idx_ == cm.null_idx ) = [];                  % never remove ( tau, p ) == ( 0, 0 )
            
            m_occ = sum( cm.m( :, occ_idx_ ) .* conj( cm.m( :, occ_idx_ ) ) );
            b_del = m_occ <= cm.epsilon^2;
            n_del = sum( b_del );
            
            if ( cm.verbose )
                
                fprintf( 1, '\t-\t %d', sum( b_del ) );
                
            end
            
            if ( n_del > 0 )
                
                cm.m( :, occ_idx_( b_del ) ) = 0;
                cm.occ( occ_idx_( b_del ) ) = false;
                cm.n_occ = cm.n_occ - n_del;
                cm.occ_idx = find( cm.occ );
                
            end
            
            if ( cm.verbose )
                
                fprintf( 1, '\t=\t %d\n', cm.n_occ );
                
            end
            
            % update logging information
            
            cm.log.t = [ cm.log.t, cm.log.t( end ) + tau_ ];
            cm.log.n_occ = [ cm.log.n_occ, cm.n_occ ];
            cm.log.n_del = [ cm.log.n_del, n_del ];
            
        end
        
        function spoiler ( cm )
            % prepare everything, if not done yet
            
            if ( ~cm.b_prep )
                
                cm.prep;
                
            end
            
            % set all transverse magnetization to zero (instantaneously)
           
            cm.m( 1 : 3 : end, : ) = 0;
            cm.m( 3 : 3 : end, : ) = 0;
                        
        end
        
        %% get results
        
        function res = sum ( cm, param )            
            % return the CM sum
            %
            % IN
            %
            % structure param (if empty, complete sum over all stored
            % configurations is returned)
            %
            % param.occ : over which configurations shall the sum be taken?
            % (optional)
            % (if field is missing, the sum is over all stored
            % configurations)
            %
            % param.omega : bulk off-resonance [ rad / ms ] (optional)
            % (default == 0, if field is missing)
            % may be an array
            %
            % param.x : position vector [ um ] (optional)
            % (== zero vector, if field is missing)
            % may have a second dimension
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
            
            if ( ~isempty( param ) && isfield( param, 'occ' ) )
                
                occ_ = param.occ;
                n_occ_ = sum( occ_ );
                
            else
                
                occ_ = cm.occ;
                n_occ_ = cm.n_occ;
                
            end
                           
            % initialize weigthing factor
            % could also be used for external weighting (currently not
            % implemented)
            % final expected size (if not equal to one) 
            % == [ 1, n_occ_, n_om, n_x ]
            
            w_n_ = 1;

            % off-resonance(s)
            
            if ( ~isempty( param ) && isfield( param, 'omega' ) )
                
                n_om = length( param.omega ); 
                omega_ = reshape( param.omega, [ 1, 1, n_om ] );
                
                % update weights
                % size( w_n_ ) == [ 1, n_occ_, n_om ]
                
                w_n_ = reshape( exp( ( - 1i .* reshape( cm.tau( occ_ ), [ n_occ_, 1 ] ) ) .* ...
                    reshape( omega_, [ 1, n_om ] ) ), [ 1, n_occ_, n_om ] );
            
            end

            % location(s)
            
            if ( ~isempty( param ) && isfield( param, 'x' ) )
                
                px_ = sum( ...
                    reshape( param.x, 3, 1, 1, [] ) .* ...
                    full( cm.p( :, occ_ ) ) );
                
                w_n_ = w_n_ .* exp ( - 1i .* px_ );
                
            end
            
            % effects due to inhomogeneous broadening (if applicable, most commonly due to R2p)
            
            if ( ~isempty( cm.inhomogeneous_decay ) )
                
                w_n_ = w_n_ .* cm.inhomogeneous_decay( cm.tau( occ_ ), param );
                
            end
            
            % calculate the isochromat
            
            res.xy = CoMo.sqrt_2 .* ...
                sum( w_n_ .* cm.m( 1 : 3 : end, occ_ ), [ 1, 2 ] );
            
            res.z = ...
                sum( w_n_ .* cm.m( 2 : 3 : end, occ_ ), [ 1, 2 ] );
                        
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
            
            if ( cm.d_tau == 0 )
                
                error( 'd_tau must be positive.' );
                
            end                        
            
            % check proton densities
            
            if ( size( cm.pd, 2 ) ~= cm.n_tissues )
                
                error( 'Proton density not set properly.' );
                
            end
            
            % check correct long time limit of magnetization transfer/exchange
            % (particle conservation is properly incorporated by set.k) 
            
            if ( ~isempty( cm.k ) )
               
                % cm.pd must lie in the null space of cm.k:
                
                k_pd = cm.k * cm.pd( : );
                
                if ( max( abs( k_pd ) > 10 * eps ) )
                    
                    % Maybe there is an error, maybe it is just numerical
                    % inaccuracy. Provide an informal message to the user and let
                    % the program continue.
                    
                    fprintf( 1, 'Check, whether cm.pd lies in null space of cm.k\n' );
                    fprintf( 1, 'Found result for cm.k * cm.pd =\n' );
                    disp( k_pd );
                    
                end
                
            end
            
        end
        
        function prep ( cm )
            % Initializes magnetization and other stuff
            
            % check, whether required parameters have been set properly
            
            cm.check;
            
            % no diffusion AND magnetization exchange
            
            if ( cm.diffusion ~= CoMo.No_Diff && ~isempty( cm.k ) )
                
                error( 'Combination of diffusion and magnetization exchange not supported.')
                
            end
            
            %% initialize everything
            
            % occupation of configuratin space
            
            cm.n_occ = 1;                    
            cm.occ_idx = 1;
            cm.null_idx = 1;
            
            cm.occ = false( 1, cm.alloc );
            cm.occ( cm.null_idx ) = true;                                  
            cm.occ_sub = zeros( cm.alloc, 4 );

            % dimensions of configuration space
            
            cm.tau = zeros( 1, cm.alloc );
            cm.p = zeros( 3, cm.alloc );

            % magnetization density is assumed to be in thermal
            % equilibrium, initially
            
            cm.m_dim = 3 * cm.n_tissues;
            cm.m = zeros( cm.m_dim, cm.alloc );
            cm.m( 2 : 3 : end, cm.null_idx ) = cm.pd;
            
            % prepare logging information
            
            cm.log.t = 0;
            cm.log.n_occ = 0;
            cm.log.n_del = 0;
            
            % everything is prepared now
            
            cm.b_prep = true;
            
        end   
        
    end
    
end