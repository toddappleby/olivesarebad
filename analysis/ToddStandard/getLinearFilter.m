function [lfilter,mfilter,pp_exp,negLexp,Cexp] = getLinearFilter(stimulus, response, varargin)

    ip = inputParser();
    ip.addParameter('analysisType', 'mle', @(x)ischar(x));
    ip.addParameter('fourierCorrection', false, @(x)islogical(x));
    ip.addParameter('binRate',1e3,@(x)isfloat(x));
    ip.addParameter('filterTime',0.5,@(x)isfloat(x));
    ip.addParameter('frameRate',60.0,@(x)isfloat(x));
    ip.parse(varargin{:});

    binRate = ip.Results.binRate;
    analysisType = ip.Results.analysisType;
    fourierCorrection = ip.Results.fourierCorrection;
    filterTime = ip.Results.filterTime;
    frameRate = ip.Results.frameRate;

    % Check inputs.
    if size(response,1) > size(response,2)
        response = response';
    end

    if size(stimulus,1) > size(stimulus,2)
        stimulus = stimulus';
    end

    % Convert filter time to sample points.
    filterPts = round(filterTime*binRate);

    if fourierCorrection
        lfilter = getLinearFilterOnline([zeros(size(response,1),round(0.1*binRate)) stimulus zeros(size(response,1),round(0.1*binRate))],...
            [zeros(size(response,1),round(0.1*binRate)) response zeros(size(response,1),round(0.1*binRate))],binRate,floor(frameRate/2));
    else
        lfilter = mean((fft([zeros(size(response,1),round(0.1*binRate)) response zeros(size(response,1),round(0.1*binRate))],[],2).*conj(fft([zeros(size(response,1),round(0.1*binRate)) stimulus zeros(size(response,1),round(0.1*binRate))],[],2))),1);
        lfilter = real(ifft(lfilter));
    end

    % Make it a unit vector.
    disp(norm(lfilter))
    lfilter = lfilter / norm(lfilter);

    if strcmpi(analysisType,'mle')
        % Need to transpose stimulus and response before sending it along.
        stimulus = stimulus';
        response = response';
        [mfilter,pp_exp,negLexp,Cexp] = getMLE(stimulus, response/binRate, flipud(lfilter(1:filterPts)'), binRate, filterPts);
    elseif strcmpi(analysisType,'istac')
        % Do iSTAC fitting of linear filters.
        stimulus = stimulus';
        response = response';
        
        [sta, stc, rawmu, rawcov] = simpleSTC(stimulus(:),response(:)-min(response(:)),filterPts);  % compute STA and STC

        % Compute iSTAC estimator
        nistacFilts = 3; % number of iSTAC filters to compute
        [mfilter,vals,DD] = compiSTAC(sta(:),stc,rawmu,rawcov,nistacFilts); % find iSTAC filters

        % Fit iSTAC model nonlinearity using filters 1 and 2
%         pp_istac12 = fitNlin_expquad_ML(stimulus,response,istacFilts(:,1:2),binRate); % LNP model struct
        pp_exp = [];
        negLexp = [];
        Cexp = [];
    else
        mfilter = [];
        pp_exp = [];
        negLexp = [];
        Cexp = [];
    end
end

% Maximum likelihood estimation of linear filter.
function [mfilter,pp_exp,negLexp,Cexp] = getMLE(stimulus, response, lfilter, binRate, filterPts)
    pp0 = makeFittingStruct_LNP(lfilter(:),binRate,[]); 

    ktbasprs.neye = 0; % number of "identity"-like basis vectors
    ktbasprs.ncos = 16; %8; % number of raised cosine basis vectors
    ktbasprs.kpeaks = [0 filterPts/2+3]; % location of 1st and last basis vector bump [0 filterPts/2+3];
    ktbasprs.b = 7; % determines how nonlinearly to stretch basis (higher => more linear)
    [ktbas, ktbasis] = makeBasis_StimKernel(ktbasprs, filterPts); % make basis
    filtprs_basis = (ktbas'*ktbas)\(ktbas'*lfilter);  % filter represented in new basis
    sta_basis = ktbas*filtprs_basis;

    % Insert filter basis into fitting struct
    pp0.k = sta_basis; % insert sta filter
    pp0.kt = filtprs_basis; % filter coefficients (in temporal basis)
    pp0.ktbas = ktbas; % temporal basis
    pp0.ktbasprs = ktbasprs;  % parameters that define the temporal basis

    % negL0 = -logli_LNP(pp0,stim(:),resp(:)); 

    opts = {'display', 'off', 'maxiter', 1000};
    [pp_exp,negLexp,Cexp] = fitLNP_1filt_ML(pp0,stimulus(:),response(:),opts);
    mfilter = pp_exp.k / norm(pp_exp.k);
end
