function [xi_phiMin,phi_phiMin,xi_phiMax,phi_phiMax,xi_xiMax,phi_xiMax] = ComputePowerControlLimits(Pdes,h,etas,phiSpec,PhiSpec,p,s,f)
% This function only looks at phi/Phi limits during nominal operation, when
% phi and Phi belong to the interval (0,1). To be expanded in future
% revisions.
%% Input 

% Sampling settings
N = s.Env.Nphi;
M = s.Env.Nxi;
phi_array = linspace(1e-9,1-1e-9,N);
xi_array = linspace(1e-9,1-1e-9,M);

% Compute available power in current flight condition
if strcmp(p.config,'e-1') || strcmp(p.config,'dual-e')
    P_avail = Pdes.EM1;
elseif strcmp(p.config,'e-2')
    P_avail = Pdes.EM2;
else
    P_avail = Pdes.GTM*f.Alpha(f.rho(h));
end

% Powertrain configuration, determines which variables have to be analyzed
switch p.config
    case 'conventional'; AnalysisCase = 0;
    case 'turboelectric'; AnalysisCase = 0;
    case 'serial'; AnalysisCase = 1;
    case 'parallel'; AnalysisCase = 1;
    case 'PTE'; AnalysisCase = 2;
    case 'SPPH'; AnalysisCase = 3;
    case 'e-1'; AnalysisCase = 0;
    case 'e-2'; AnalysisCase = 0;
    case 'dual-e'; AnalysisCase = 2;
end

% For the SPPH, either a given phi value or a given Phi value have to be
% specified (not both)
if AnalysisCase == 3
    
    % Set input phiSpec and PhiSpec to NaN if unspecified
    if isempty(phiSpec); phiSpec = NaN; end
    if isempty(PhiSpec); PhiSpec = NaN; end
    
    % Check input and assign analysis method
    % If none of the two have been specified
    if isnan(phiSpec) && isnan(PhiSpec)
        error(['A value has to be assigned to either phi or Phi, '...
                'at which the envelope is evaluated '...
                'for a SPPH configuration.'])
            
    % If both have been specified
    elseif ~isnan(phiSpec) && ~isnan(PhiSpec)
        error(['Only one value can be assigned (either phi or Phi) '...
                'at which the envelope is evaluated '...
                'for a SPPH configuration.'])        
    end
end


%% Compute envelopes

% For configurations where only throttle is used
if AnalysisCase == 0
    
    % phi and Phi limits can be ignored
    xi_phiMin = NaN;
    phi_phiMin = NaN;
    xi_phiMax = NaN;
    phi_phiMax = NaN;
    phi_xiMax = NaN;
    xi_xiMax = NaN;
    
% Compute envelope for powertrains with at least two DOFs 
% Note: this segment can be evaluated for xi and phi, or for xi and Phi. To
% reduce script length, Phi is also named "phi". But it refers to shaft
% power ratio nonetheless! (for PTE or dual-e powertrains)
else    
    
    % Initialize arrays
    phi_phiMin = NaN(M,1);
    phi_phiMax = NaN(M,1);
    xi_phiMin = NaN(M,1);
    xi_phiMax = NaN(M,1);
    toggleGT = zeros(M,N); % Toggle arrays are useful for debugging
    toggleEM1 = zeros(M,N); 
    toggleEM2 = zeros(M,N); 
    togglebat = zeros(M,N); 
    toggles1 = zeros(M,N);
    P_out = cell(M,N);
        
    % Evaluate powertrain M x N times to get response at each phi and xi.
    % For each xi
    for i = 1:M
        
        % Reset variables
        phi_out = NaN(1,N);
        
        % Loop over phi values to get limit for the given xi
        for j = 1:N

            % For powertrains with only phi
            if AnalysisCase == 1
                [P_out{i,j},~,~,~,~] = ...
                   PowerTransmissionComputation_v2(p.config,...
                   etas,phi_array(j),[],xi_array(i),[],[],[],P_avail);
                
            % For powertrains with only Phi
            elseif AnalysisCase == 2;
                [P_out{i,j},~,~,~,~] = ...
                   PowerTransmissionComputation_v2(p.config,...
                   etas,[],phi_array(j),xi_array(i),[],[],[],P_avail);
                
            % For powertrains with both phi and Phi, when a specific phi
            % value has been specified
            elseif AnalysisCase == 3 && ~isnan(phiSpec);
                [P_out{i,j},~,~,~,~] = ...
                   PowerTransmissionComputation_v2(p.config,...
                   etas,phiSpec,phi_array(j),xi_array(i),[],[],[],P_avail);
                
            % For powertrains with both phi and Phi, when a specific Phi
            % value has been specified
            elseif AnalysisCase == 3 && ~isnan(PhiSpec);
                [P_out{i,j},~,~,~,~] = ...
                   PowerTransmissionComputation_v2(p.config,...
                   etas,phi_array(j),PhiSpec,xi_array(i),[],[],[],P_avail);
            end
                        
            % Manually select the power paths which we want to compare to
            % powers obtained from the WPdes array in the power loading
            % diagram. Select one path for each of the five branches in the
            % powertrain architecture model. If on a given branch one path is
            % above the limit, then any other path on that branch will also be.
            PGT = P_out{i,j}.gt;
            PEM1 = max([abs(P_out{i,j}.gb) abs(P_out{i,j}.e1)]);
            PEM2 = max([abs(P_out{i,j}.e2) abs(P_out{i,j}.s2)]);
            Pbat = abs(P_out{i,j}.bat);
            Ps1 = abs(P_out{i,j}.s1);
            
            % Check if they surpass the design power
            if PGT > Pdes.GTM; toggleGT(i,j) = 1; end
            if PEM1 > Pdes.EM1; toggleEM1(i,j) = 1; end
            if PEM2 > Pdes.EM2; toggleEM2(i,j) = 1; end
            if Pbat > abs(Pdes.bat); togglebat(i,j) = 1; end
            if Ps1 > abs(Pdes.s1); toggles1(i,j) = 1; end
            
            % Gather toggles
            indexes = sum([toggleGT(i,j) toggleEM1(i,j) toggleEM2(i,j)...
                           togglebat(i,j) toggles1(i,j)]);
            
            % Save current phi value if limit is not surpassed
            if indexes == 0
                phi_out(j) = phi_array(j);
            else
                phi_out(j) = NaN;
            end
            
            % Break loop if minimum and maximum have both been found to
            % save time
            if j > 2
                if ~isnan(phi_out(1)) && isnan(phi_out(j))
                    break
                elseif isnan(phi_out(1)) && isnan(phi_out(j)) && ...
                        ~isnan(phi_out(j-1))
                    break
                end
            end
        end
        
        % Select limiting phi values
        phi_phiMax(i) = max(phi_out);
        xi_phiMax(i) = xi_array(i);
        phi_phiMin(i) = min(phi_out);
        xi_phiMin(i) = xi_array(i);
    end
    
    % Add line for possible GT limitation
    if any(isnan(phi_phiMax))
        [xi_xiMax,idx] = max(xi_array(~isnan(phi_phiMax)));
        xi_xiMax = [xi_xiMax xi_xiMax];
        phi_xiMax = [phi_phiMin(idx) phi_phiMax(idx)];
    else
        xi_xiMax = [NaN NaN];
        phi_xiMax = [NaN NaN];
    end
end






