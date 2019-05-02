function pmf_FastOddball_Numerosity( varargin )
    addpath(genpath('./common')); % add xDiva common folder

    % _VV_2015_0206 so that 'definitions' are always initialized
    definitions = MakeDefinitions;
    parameters	= {};			% always varargin{2}			initialize here for sceop
    timing		= {};			% always varargin{3}
    videoMode	= {};			% always varargin{4}
    iS = 1;		% indices of parts S,B,1,2 in definitions cell array
    iB = 2;
    i1 = 3;
    i2 = 4;
    iSweep = 5;
    iMod   = 6;
    % 	iAux   = 7;

    if nargin == 0
        TestSubFunction			% does xDiva ever call w/o inputs or was this just for debugging?
        % 		error('No input')
    elseif ismember( varargin{1}, { 'GetDefinitions', 'ValidateParameters', 'MakeMovie', 'ValidateDefinition', 'TestSubFunction' } )
        eval( varargin{1} )		% feval throws 'Undefined function or variable' error.  need to suppress output?
    end
    
    return

    function rV = MakeDefinitions
        % - currently implementing 'integer', 'double', 'nominal'
        % - types of the cells in each parameter row
        % - only "standard" type names can be used in the "type" fields
        % - 'nominal' params should have
        %     (a) at least one item
        %		(b) value within the size of the array
        % - all other params (so far) should have empty arrays of items
        
        rV = { ...
            
        
        % - Part 'S' parameters must use standard parameter names
        % - 'Sweep Type' : at least 'Fixed' sweep type must be defined
        % - 'Modulation' : at least 'None' modulation type must be defined
        % - 'Step Type'  : at least 'Lin Stair' type must be defined,
        %                  first 4 step types are reserved, custom step types can only be added after them
        {
        'View Dist (cm)'	57.0		'double'		{}
        'Mean Lum (cd)'		50.0		'double'		{}
        %			'Fix Point'			 1.0		'double'		{}	% 2016.07.22 changed to 'nominal'
        'Fix Point'			'None'		'nominal'	{ 'None', 'Cross', 'Up Cue', 'Down Cue', 'Left Cue', 'Right Cue', 'Nonius', 'Block A', 'Block B', 'Fancy' }
        'Sweep Type'		'Fixed'		'nominal'	{ 'Fixed', 'Contrast' }
        'Step Type'			'Lin Stair'	'nominal'	{ 'Lin Stair', 'Log Stair', 'Lin Sin', 'Log Sin' }
        'Sweep Start'		 0.0		'double'		{}
        'Sweep End'			 1.0		'double'		{}
        'Modulation'		'None'		'nominal'	{ 'None', 'Square', 'Sine' }
        }
        
        
        % - Part 'B' definition cell must be defined (at least as an empty array)
        {
        'ModInfo'			1.0				'double'		{}
        'Target Gamma'		1.8				'double'		{}
        'Numerosity covaries with ' 'Random' 'nominal'	{ 'Area', 'Size', 'Density', 'Random'}
        'Control Range: Lowest'     '5'          'nominal'	{ '1','2','3','4','5','6','7','8','9' }
        'Control Range: Highest'	'9'          'nominal'	{ '1','2','3','4','5','6','7','8','9' }
        'Control Type'              'Full'        'nominal'	{ 'Full','Limited' }
        'Stimulus Size (dva x dva)'	    10.0	  'double'		{}
        %            'Use Fill Method'   'No'            'nominal'   { 'Yes', 'No'}
        %			'Scale'				'1x'			'nominal'	{ '1x', '2x' }				% VV's going to build this into xDiva
        % 			'V Size (deg)'		2.0				'double'		{}
        }
        
        
        % - Part '1' 'notused' allows to skip creating GUI element in particular column
        %            'Cycle Frames' has to be 1st row
        
        {
        'Cycle Frames'		60.0			'integer'	{}
        % 			'Image Category'	'Faces'		'nominal'	{ 'Faces', 'Objects' }
        'Contrast (rel)'	1.0			'double'		{}
        'Randomization'	'Yes'			'nominal'	{ 'Yes', 'No' }
        'Operation'			'None'		'nominal'	{ 'None', 'H Flip', 'V Flip', 'Invert' }
        'Numerosity'         '5'        'nominal'	{ '1','2','3','4','5','6','7','8','9' }
        'Shape'              'circle'   'nominal'	{ 'circle','gabor'} % 'triangle','square'
        }
        
        % - Part '2' ?optional?
        {
        'Cycle Frames'		12.0		'integer'	{}
        % 			'Image Category'	'Objects'	'nominal'	{ 'Faces', 'Objects' }
        'Contrast (rel)'	1.0			'double'	{}
        % 			'Contrast (pct)'	90.0		'double'	{}
        'Randomization'		'Yes'		'nominal'	{ 'Yes', 'No' }
        'Operation'			'None'		'nominal'	{ 'None', 'H Flip', 'V Flip', 'Invert' }
        'Numerosity'         '6'        'nominal'	{ '1','2','3','4','5','6','7','8','9' }
        'Shape'              'circle'   'nominal'	{ 'circle','gabor'} % 'triangle','square'
        }
        
        
        
        
        % Sweepable parameters
        % The cell array must contain as many rows as there are supported Sweep Types
        % 1st column (Sweep Types) contains Sweep Type as string
        % 2nd column (Stimulus Visiblity) contains one of the following strings,
        % indicating how stimulus visibility changes when corresponding swept parameter value increases:
        %   'constant'   - stimulus visibility stays constant
        %   'increasing' - stimulus visibility increases
        %   'decreasing' - stimulus visibility decreases
        % 3rd column contains a single-row cell array of pairs, where each pair is a single-row cell
        % array of 2 strings: { Part name, Parameter name }
        % If sweep affects only one part, then you only need one { part, param} pair.
        % If it affects both parts, then you need both  pairs, e.g. for "Contrast" below
        % Standard part names: 'S', 'B', '1', '2'
        {
        'Fixed'			'constant'		{}
        'Contrast'		'increasing'	{ { '1', 'Contrast (rel)' }, { '2', 'Contrast (rel)' } }		% everything in list gets swept
        }
        
        % ModInfo information
        % The cell array must contain as many rows as there are supported Modulations
        % 1st column (Modulation) contains one of the supported Modulation typs as string
        % 2nd column contains the name of the ModInfo parameter as string
        % 3rd column (default value) contains default value of the ModInfo
        % parameter for this Modulation
        {
        'None'				'ModInfo'			 0.0
        'Square'			'OffLum (Cd)'		50.0
        'Sine'				'OffLum (Cd)'		50.0
        }
        
        % Auxiliary information.
        % Required by xDiva, but not by Matlab Function
        {
        'Version'						1				% of the matlab function?
        'Adjustable'					true
        'Needs Unique Stimuli'		true			% e.g. if you want to re-randomize something every trial set to true
        'Supports Interleaving'		true
        'Part Name'						{ 'Oddball'  'Reference' }
        'Frame Rate Divisor'		{ 4 2 }
        'Max Cycle Frames'			{ 120 30 }
        'Allow Static Part'			{ false false }
        % 			'Supported Psy Procedures'	{}				% currently not used
        % 			'Supported Psy Nulls'		{}				% currently not used
        }
        
        };
    
    end

    function GetDefinitions
        assignin( 'base', 'output', definitions );
    end

    function ValidateParameters
        % xDiva invokes Matlab Engine command:
        % pmf_<subParadigmName>( 'ValidateParameters', parameters, timing, videoMode );
        % "parameters" here is an input argument. Its cellarray has the same structure
        % as "definitions" but each parameter row has only first two elements
        
        % The "timing" and "videoMode" cellarrays have the same row
        % structure with each row having a "name" and "value" elements.
        
        % timing names from TimingRec declaration in PowerDiva.h
        % videoMode names from VideoModeRx declaration in VxVideoMode.h
        
        
        % nFrameCycle1 = 1x or more multiple of 2 if square wave
        %              = 1x or more multiple of 4 if sine wave
        % nFrameCycle2 = 2x or more multiple of nFrameCycle1
        % nFrameStep / nFrameCycle2 = integer >= 1
        % nFrameBin  / nFrameCycle2 = integer >= 1
        
        [ parameters, timing, videoMode ] = deal( varargin{2:4} );
        % save parameters to desktop, for debugging
        save('~/Desktop/xDiva_params.mat','parameters', 'timing','videoMode');
        
        nFrameCycle1 = xDiva.GrabCellValue( parameters{i1}, 'Cycle Frames' );			% double
        nFrameCycle2 = xDiva.GrabCellValue( parameters{i2}, 'Cycle Frames' );
        nFrameStep	 = xDiva.GrabCellValue( timing, 'nmbFramesPerStep' );
        nFrameBin	 = xDiva.GrabCellValue( timing, 'nmbFramesPerBin' );
        
        validationMessages = {
            'frames/cycle1 not multiple of 2'
            'frames/cycle1 not multiple of 4'
            'frames/cycle2 not multiple of frames/cycle1'
            'frames/step not multiple of frames/cycle2'
            'frames/bin not multiple of frames/cycle2'
            };
        parValidFlags = true(size(validationMessages));
        switch xDiva.GrabCellValue( parameters{iS}, 'Modulation' )
            case 'None'
            case 'Square'
                parValidFlags(1) = mod(nFrameCycle2,2) == 0;
            case 'Sine'
                parValidFlags(2) = mod(nFrameCycle2,4) == 0;
        end
        parValidFlags(3) = nFrameCycle1 >= 2*nFrameCycle2 && mod( nFrameCycle1, nFrameCycle2 ) == 0;
        parValidFlags(4) = nFrameStep   >=   nFrameCycle1 && mod( nFrameStep  , nFrameCycle1 ) == 0;
        parValidFlags(5) = nFrameBin    >=   nFrameCycle1 && mod( nFrameBin   , nFrameCycle1 ) == 0;
        parValidFlag = all( parValidFlags );
        if ~parValidFlag
            validationMessages = validationMessages(~parValidFlags);
        else
            validationMessages = cell(1);
        end
        
        % validate and correct range options
        numero(1) = str2double(xDiva.GrabCellValue( parameters{i1}, 'Numerosity' ));
        numero(2) = str2double(xDiva.GrabCellValue( parameters{i2}, 'Numerosity' ));
        mindesired = str2double(xDiva.GrabCellValue( parameters{iB}, 'Control Range: Lowest' ));
        maxdesired = str2double(xDiva.GrabCellValue( parameters{iB}, 'Control Range: Highest' ));
        if mindesired >= maxdesired
            controlStr = sprintf('Control Range Lowest %d must be lower than Control Range Highest %d.',mindesired,maxdesired);
            count = 0;
            while mindesired >= maxdesired
                count = count +1;
                if mindesired > 1 && mod(count,2)
                    mindesired = mindesired - 1;
                elseif maxdesired < 9 && ~mod(count,2)
                    maxdesired = maxdesired + 1;
                else
                end
            end
            controlStr = sprintf('%s Control Range adjusted to %d-%d.',controlStr,mindesired,maxdesired);
            validationMessages = xDiva.AppendVMs(validationMessages,controlStr);
            parValidFlag = false;
        else
        end
        if sum(ismember(numero,mindesired:maxdesired)) ~= 2;
            controlStr = 'Control Range must include both reference and oddball.';
            minDiff = mindesired-numero;
            maxDiff = maxdesired-numero;
            if sum(minDiff > 0); % if any numero is smaller than lowest range
                controlStr = sprintf('%s Control Range Lowest %d adjusted to %d.',controlStr,mindesired,min(numero));
                mindesired = min(numero); % assign lowest numero
            else
            end
            if sum(maxDiff < 0); % if any numero is bigger than highest range
                controlStr = sprintf('%s Control Range Highest %d adjusted to %d.',controlStr,maxdesired,max(numero));
                maxdesired = max(numero);
            else
            end
            validationMessages = xDiva.AppendVMs(validationMessages,controlStr);
            parValidFlag = false;
        else
        end
        parameters = xDiva.ReplaceParam(parameters, 'B', 'Control Range: Lowest', num2str(mindesired));
        parameters = xDiva.ReplaceParam(parameters, 'B', 'Control Range: Highest', num2str(maxdesired));
        
        if numero(1)/numero(2) < 0.5 && strcmp(xDiva.GrabCellValue( parameters{iB}, 'Control Type' ),'Full')
            controlStr = sprintf('Oddball %d is less than half of Reference %d.\n',numero(1),numero(2));
            controlStr = sprintf('%sControl Type cannot be "Full" - set to "Limited".\n',controlStr);
            controlStr = sprintf('%sNote that more complete control of non-number parameters can be achieved with different number sets\n',controlStr);
            validationMessages = xDiva.AppendVMs(validationMessages,controlStr);
            parValidFlag = false;
            parameters = xDiva.ReplaceParam(parameters, 'B', 'Control Type', 'Limited');
        else
        end
        
        % check that stimulus size is not bigger than screen
        stimSize = xDiva.GrabCellValue( parameters{iB}, 'Stimulus Size (dva x dva)' ) * 60; % stim size in arc minutes
        [pix2am,maxWidthAM,maxHeightAM] = xDiva.Pix2Arcmin(videoMode,parameters);
        
        if stimSize > min(maxWidthAM,maxHeightAM);
            [newStimSize,stimIdx] = min([maxWidthAM,maxHeightAM]);
            if stimIdx == 1
                controlStr = sprintf('Stimulus size %d exceeds the display width.\n',stimSize/60);
            else
                controlStr = sprintf('Stimulus size %d exceeds the display height.\n',stimSize/60);
            end
            controlStr = sprintf('%sStimulus size adjusted to maximum possible, %.2f.\n',controlStr,newStimSize / 60);
            validationMessages = xDiva.AppendVMs(validationMessages,controlStr);
            parValidFlag = false;
            parameters = xDiva.ReplaceParam(parameters, 'B', 'Stimulus Size (dva x dva)' , newStimSize / 60);
        else
        end
        
        % validate and correct gabor options
        %if strcmp(xDiva.GrabCellValue( parameters{iB}, 'Use Fill Method' ),'Yes')
        %    gaborOpts = [strcmp(xDiva.GrabCellValue( parameters{i1}, 'Shape' ),'gabor'),strcmp(xDiva.GrabCellValue( parameters{i2}, 'Shape' ),'gabor')];
        %    if sum(gaborOpts)
        %        validationMessages = xDiva.AppendVMs(validationMessages,'Fill Method cannot be used with Gabors, setting "Use Fill Method" to "No"');
        %        parameters = xDiva.ReplaceParam(parameters, 'B', 'Use Fill Method', 'No');
        %        parValidFlag = false;
        %    else
        %    end
        %else
        %end
        assignin( 'base', 'output', { parValidFlag, parameters, validationMessages } )
    end

    function ValidateDefinition
        % Test to make sure that 'definitions' array is correct, i.e.
        % correct number of parts, structure of param items, etc. xDiva
        % will call this first to ensure that parameters contain no
        % unhandled errors.  This can also be used by Matlab user during
        % paradigm development.
        
        % Note, the 'definitions' array already exists
        
        % ... check definitions here ...
        tIsValidDefinition = true;
        
        assignin( 'base', 'validDefinition', tIsValidDefinition );
    end

    function MakeMovie
        % xDiva creates variables "parameters", "timing", "videoMode" in
        % Matlab Engine workspace and invokes Matlab Engine command:
        
        % pmf_<subParadigmName>( 'MakeMovie', parameters, timing, videoMode );
        % e.g.,
        % pmf_MyBilateralGrating( 'MakeMovie', parameters, timing, videoMode );
        
        % where pmf_<subParadigmName> is the name of the m-file selected by
        % the MatlabFunction paradigm dialog "Choose" control
        
        % Don't bother with fixation point in actual movie, xDiva will add
        % it later.
        
        %		*** Loading precalculated images from file ***
        %		tPFNm = '/xDiva/xDiva_Data/MatlabFiles_Gliders_Clouds/glider/pilotStimGliders20Right_F_001.mat';
        %		tD = load( tPFNm, 'images', 'imageSequence' );
        % 		assignin( 'base', 'images', tD.images );
        % 		assignin( 'base', 'imageSequence', tD.imageSequence );
        
        %try
            [ parameters, timing, videoMode, trialNumber ] = deal( varargin{2:5} );
            
            nFrameBin     = xDiva.GrabCellValue( timing, 'nmbFramesPerBin' );
            nFrameStep    = xDiva.GrabCellValue( timing, 'nmbFramesPerStep' );
            nBinPrelude   = xDiva.GrabCellValue( timing, 'nmbPreludeBins' );
            nStepCore     = xDiva.GrabCellValue( timing, 'nmbCoreSteps' );
            nFramePrelude = nFrameBin * nBinPrelude;
            nFrameCore    = nFrameStep * nStepCore;
            nFrameTrial   = nFrameCore + 2 * nFramePrelude;
            
            modType = xDiva.GrabCellValue( parameters{iS}, 'Modulation' );
            if ~strcmp( modType, 'None' )
                lumOff = xDiva.GrabCellValue( definitions{iMod}, modType, 3 );
            else
            end
            
            nFrameCycleOdd = xDiva.GrabCellValue( parameters{i1}, 'Cycle Frames' );       % double
            nFrameCycleRef = xDiva.GrabCellValue( parameters{i2}, 'Cycle Frames' );
            frameRatioOddRef = nFrameCycleOdd / nFrameCycleRef;
            
            preludeType = xDiva.GrabCellValue( timing, 'preludeType' );                   % 0=dynamic, 1=blank, 2=static
            nImgOdd = nFrameCore / nFrameCycleOdd;									% # images to show per trial
            nImgRef = nImgOdd * ( frameRatioOddRef - 1 );
            
            if preludeType == 0		% dynamic
                nCycle2Prelude = 2 * nFramePrelude / nFrameCycleOdd;
                nImgOdd = nImgOdd + nCycle2Prelude;
                nImgRef = nImgRef + nCycle2Prelude * ( frameRatioOddRef - 1 );
            else
            end
                        
            %% THIS IS  THE NUMEROSITY CODE BEGINS
                        
            stimSize = xDiva.GrabCellValue( parameters{iB}, 'Stimulus Size (dva x dva)' ) * 60; % stim size in arc minutes
            windowSize = floor(stimSize / xDiva.Pix2Arcmin(videoMode,parameters)); % stim size in pixels
            
            numeroOdd = str2double(xDiva.GrabCellValue( parameters{i1}, 'Numerosity' ));
            numeroRef = str2double(xDiva.GrabCellValue( parameters{i2}, 'Numerosity' ));
            shapeOdd = xDiva.GrabCellValue( parameters{i1}, 'Shape' );
            shapeRef = xDiva.GrabCellValue( parameters{i2}, 'Shape' );
            covary = find(cell2mat(cellfun(@(x) strcmp(x,xDiva.GrabCellValue( parameters{iB}, 'Numerosity covaries with ' )),{ 'Area', 'Size', 'Density', 'Random'},'uni',false)));
            minDesired = str2double(xDiva.GrabCellValue( parameters{iB}, 'Control Range: Lowest' ));
            maxDesired = str2double(xDiva.GrabCellValue( parameters{iB}, 'Control Range: Highest' ));
            fullControl = strcmp(xDiva.GrabCellValue( parameters{iB}, 'Control Type' ),'Full');
            useFill = false; %strcmp(xDiva.GrabCellValue( parameters{iB}, 'Use Fill Method' ),'Yes');
            
            % initialize 4D .mat
            img = zeros(windowSize,windowSize,1,nImgRef+nImgOdd,'uint8');
            
            criticalItemSize = minDesired/maxDesired; % this is the critical item size defined by red arrows in accompanying explanatory graphs
            criticalToa = minDesired/maxDesired;
            rMax = min( (1/sqrt(maxDesired))/2.5 , .3); % crazy code from Dehaene to compute maximum radius of dots
            
            % generate luminance values
            lumRange = linspace(xDiva.GrabCellValue( videoMode,  'minLuminanceCd' ),xDiva.GrabCellValue( videoMode,  'maxLuminanceCd' ),256);
            lumIdx = findClosest(xDiva.GrabCellValue( parameters{iS},  'Mean Lum (cd)' ),lumRange);
            lumIdx = lumIdx - 1; % lowest values is zero, biggest is 255;
            
            %% SET UP STIMULUS PARAMETERS
            % generate random order of controls if needed
            
            % conditions == 1: sizeCovary = 0; areaCovary = 1
            % conditions == 2: sizeCovary = 1; areaCovary = 0
            % conditions == 3: sizeCovary = 0, areaCovary = 0
            
            if covary == 4 % random
                condLabels = [1,2,3];
            else
                condLabels = covary;
            end
            
            nCond = length(condLabels);
            if nCond > 1
                nReps = floor(nImgOdd/nCond);
                conditions = repmat(condLabels,1,nReps);
                conditions = [conditions condLabels(randi(nCond,1,nImgOdd-length(conditions)))]; % randomly add any additional necessary conditions
                conditions = conditions(randperm(length(conditions))); % randomize condition order
            else
                conditions = ones(1,nImgOdd)*condLabels;
            end
            
            refCond = reshape(repmat(conditions,frameRatioOddRef-1,1),1,nImgRef); % total number of reference conditions
            oddCond = [zeros(1,nImgRef), conditions]; % total number of odd conditions
            
            % now generate the stimulus parameters
            imgSizeList = zeros(1,nImgRef+nImgOdd);
            imgAreaList = zeros(1,nImgRef+nImgOdd);
            
            % variables for size and toa values, 1 & 2 = reference, 3 = oddball
            itemSize = zeros(length(condLabels),3);
            toa = zeros(length(condLabels),3);
            
            bounds = [minDesired maxDesired]; % to push parameter to a boundary when numeroOdd == numeroRef
            [~,farIdx] = max(abs(numeroRef-bounds)); %find the far condition
            for c = 1:length(condLabels)
                switch condLabels(c)
                    case 1 % area covaries, size and density does not
                        itemSize(c,1) = criticalItemSize * (minDesired/numeroRef);
                        itemSize(c,2) = criticalItemSize * (maxDesired/numeroRef);
                        itemSize(c,3) = criticalItemSize;
                        toa(c,1) = criticalToa;
                        toa(c,2) = 1;
						if numeroOdd == numeroRef
							toa(c,3) = criticalToa * (bounds(farIdx)/minDesired); % when reference is furthest from the minimum boundary, all oddball toa the same
						else
							toa(c,3) = criticalToa * (numeroOdd/minDesired); % when oddball is equal to boundary, all oddball toa same
						end
                    case 2 % size covaries, area and density does not
                        itemSize(c,1) = criticalItemSize;
                        itemSize(c,2) = 1;
                        if numeroOdd == numeroRef
                            itemSize(c,3) = criticalItemSize * (maxDesired/bounds(farIdx)); % when reference is furthest from the max boundary, all oddball size the same
                        else
                            itemSize(c,3) = criticalItemSize * (maxDesired/numeroOdd); % when oddball is equal to boundary, all oddball size same
                        end
                        toa(c,1) = criticalToa * (numeroRef/maxDesired);
                        toa(c,2) = criticalToa * (numeroRef/minDesired);
                        toa(c,3) = criticalToa;
                    case 3 % density covaries, area and size does not
                        itemSize(c,1) = criticalItemSize * (minDesired/numeroRef);
                        itemSize(c,2) = criticalItemSize * (maxDesired/numeroRef);
                        itemSize(c,3) = criticalItemSize;
                        toa(c,1) = criticalToa * (numeroRef/maxDesired);
                        toa(c,2) = criticalToa * (numeroRef/minDesired);
                        toa(c,3) = criticalToa;
                    otherwise
                end
                
                % compute the parameters for references
                refPerCond = length(find(refCond == condLabels(c)));
            
                parSqDim = ceil(sqrt(refPerCond)); % Dimension of 2-D parameter space to sample
                if ~fullControl
                    parSqDim = ceil(sqrt(refPerCond/2))*2;
                else
                end
                sizeSpace = linspace(itemSize(c,1), itemSize(c,2), parSqDim+1);
                areaSpace = linspace(toa(c,1), toa(c,2), parSqDim+1);        
                sizeRand = cell2mat(arrayfun(@(x) sizeSpace(x) + (sizeSpace(x+1) - sizeSpace(x))*rand,1:length(sizeSpace)-1,'uni',false));
                areaRand = cell2mat(arrayfun(@(x) areaSpace(x) + (areaSpace(x+1) - areaSpace(x))*rand,1:length(areaSpace)-1,'uni',false));
                [parSize, parArea] = meshgrid(sizeRand,areaRand);
                if ~fullControl
                    twoQuads = sign(parSize-(itemSize(c,1)+((itemSize(c,2)-itemSize(c,1))/2))) == sign(parArea-(toa(c,1)+((toa(c,2)-toa(c,1))/2)));
                    parSize = parSize(twoQuads);
                    parArea = parArea(twoQuads);
%                     randomreal = rand(1,refPerCond);
%                     tmpSizeList = zeros(1,refPerCond);
%                     tmpAreaList = zeros(1,refPerCond);
%                     for j = 1:length(tmpSizeList)
%                        tmpSizeList(j) = tmpSizeList(j) + itemSize(c,1) + randomreal(j)*(itemSize(c,2)-itemSize(c,1));
%                        tmpAreaList(j) = tmpAreaList(j) + toa(c,1) + randomreal(j)*(toa(c,2)-toa(c,1));
%                     end
                    
                end

                parSize = reshape(parSize,1,[]);
                parArea = reshape(parArea,1,[]);
                shflIdx = randperm(length(parSize));
                parSize = parSize(shflIdx);
                parArea = parArea(shflIdx);
                parSize = parSize(1:refPerCond);
                parArea = parArea(1:refPerCond);
                imgSizeList(refCond == condLabels(c)) = parSize;
                imgAreaList(refCond == condLabels(c)) = parArea;
%                 if ~fullControl
%                    imgSizeList(refCond == condLabels(c)) = tmpSizeList;
%                    imgAreaList(refCond == condLabels(c)) = tmpAreaList;
%                 end
                % compute the parameters for oddball
                imgSizeList(oddCond == condLabels(c)) = itemSize(c,3);
                imgAreaList(oddCond == condLabels(c)) = toa(c,3);
            end
            
            %% generate the stimuli
            global img_out; global seq_out;
            for i = 1:(nImgRef + nImgOdd)
                success = false;
                while ~success
                    if i <= nImgRef
                        [success,tempImg] = generate_set(numeroRef,shapeRef,windowSize,rMax,imgAreaList(i),imgSizeList(i),lumIdx,useFill);
                    else
                        [success,tempImg] = generate_set(numeroOdd,shapeOdd,windowSize,rMax,imgAreaList(i),imgSizeList(i),lumIdx,useFill);
                    end
                end
                img(:,:,:,i) = tempImg;
            end
            img_out = img;
            
            if length(unique(itemSize(:,3))) == 1
                fprintf('ref %d, odd %d, dot size does not vary over oddball conditions\n',numeroRef,numeroOdd);
            else
            end
            if length(unique(toa(:,3))) == 1
                fprintf('ref %d, odd %d, area does not vary over oddball conditions\n',numeroRef,numeroOdd);
            else
            end
            %% COMPUTE IMAGE SEQUENCE
            
            %tmp = load('/Users/kohler/xDiva/xDiva_ImageFiles/FastOddball/Set02GrayEq_Face.mat');
            %imgSetRef = tmp.img;
            %tmp = load('/Users/kohler/xDiva/xDiva_ImageFiles/FastOddball/Set02GrayEq_Face.mat');
            %imgSetOdd = tmp.img;
            %img = cat(4, imgSetRef, imgSetOdd );
            
            lumVals = cell2mat(arrayfun(@(x) mean(mean(img(:,:,1,x))),1:size(img,4),'uni',false));
            if max( abs(lumVals-lumIdx) ) > 2
                [~, maxIdx ] = max( abs(lumVals-lumIdx) );
                error('Luminance index (0-255) should be %0.2d, is %0.2d!',lumIdx,lumVals(maxIdx));
            else
            end
            
            clear imgSet*;
            
            % *** assuming square wave for the time being ***
            if preludeType == 0				% dynamic
                iCat1 =                1:nFrameCycleRef:nFrameTrial;
                iCat0 = nFrameCycleRef/2+1:nFrameCycleRef:nFrameTrial;
            else
                iCat1 = nFramePrelude               +1:nFrameCycleRef:nFrameTrial-nFramePrelude;
                iCat0 = nFramePrelude+nFrameCycleRef/2+1:nFrameCycleRef:nFrameTrial-nFramePrelude;
                if preludeType == 2			% static
                    iCat1(1) = 1;
                elseif preludeType == 1		% blank
                    iCat0 = [ 1, iCat0, nFrameTrial-nFramePrelude+1 ];
                else
                end
            end
            % replace onsets of img1 with img2
            iCat2 = iCat1(frameRatioOddRef:frameRatioOddRef:end);
            iCat1(frameRatioOddRef:frameRatioOddRef:end) = [];
            
            
            if ( numel(iCat1) ~= nImgRef ) || ( numel(iCat2) ~= nImgOdd )
                error('not enough images generated');
            else
            end
            
            imgSeq = zeros( [ nFrameTrial, 1 ], 'int32' );
            imgSeq(iCat1) = 1:nImgRef;
            imgSeq(iCat2) = (1:nImgOdd) + nImgRef;
            
            squareFlag = strcmp( modType, 'Square' );
            if squareFlag
                %				img(:,:,:,nImg) =  xDiva.GrabCellValue( videoMode,  'meanLuminanceBitmapValue' );		% encoded in relative luminance from min to max?
                lumMin = xDiva.GrabCellValue( videoMode, 'minLuminanceCd' );
                img(:,:,:,nImgRef + nImgOdd + 1) = (  lumOff - lumMin ) / (  xDiva.GrabCellValue( videoMode, 'maxLuminanceCd' ) - lumMin ) * 255;
                % off part of cycle
                imgSeq(iCat0) = nImgRef + nImgOdd + 1;
            else
            end
            
            assignin( 'base', 'output', { true, img, imgSeq } )			% put local var into global space as 'output'
            
            seq_out = imgSeq;
        %catch ME
        %    disp(ME.message)
        %    disp(ME.stack(1))
        %    % 			disp(ME.stack(end))
        %    assignin( 'base', 'output', { false, zeros([1 1 1 1],'uint8'), 1 } )
        %end
        
        return
        
        %{

%	--- Timing
		isBlankPrelude	= ( preludeType == 1);


%	--- Video Mode
 		wPix = xDiva.GrabCellValue( videoMode,  'widthPix' );
		hPix = xDiva.GrabCellValue( videoMode,  'heightPix' );
		wCm  = xDiva.GrabCellValue( videoMode,  'imageWidthCm' );
		hCm  = xDiva.GrabCellValue( videoMode,  'imageHeightCm' );

		% Parameter names and their string values must be exactly
		% the same as defined in MakeDefinitions subfunction
		viewDistCm	= xDiva.GrabCellValue( parameters{iS}, 'View Dist (cm)' );
		sweepType	= xDiva.GrabCellValue( parameters{iS}, 'Sweep Type' );
		isSwept		= ~strcmpi( sweepType,'Fixed' );
		sweepStart	= xDiva.GrabCellValue( parameters{iS}, 'Sweep Start' );
		sweepEnd	   = xDiva.GrabCellValue( parameters{iS}, 'Sweep End' );
        %}
        
    end


    function [success,img] = generate_set(numerosity,shape,windowSize,rmax,totaloccupiedarea,itemsize,lumIdx,useFill)
        %%%% This function is used to create a stimulus for numerosity experiments
        % S. Dehaene Version as of 12 May 2004
        %
        % This function will generate one set of objects in figure 101.
        % Arguments:
        % -	numerosity = numerosity
        % -	shape = a standardized vector defining the shape
        % -	windowSize = size of figure in pixels (the figure will be square)
        % -	rmax = maximum radius of each item (set to 0 if you want automatic
        % setting of this value (see below))
		% -	totaloccupiedarea (expressed as percentage of max, from 0 to 1)
        % -	itemsize (expressed as percentage of max, from 0 to 1)
        % - luminance index of background (min: 0, max: 255).
        % Typical example:
        %
        % generate_set(32,shape,350,0,1,0.2);
        %
        % creates a set of 32 shapes in a 350-pixel window with up to 100
        % positions for objects. The objects are distributed over the entire
        % window (1=100% of windowSize), and are at 20% of the maximum allowable
        % size
        %
        % The shape must be defined in the following way by the calling program:
        %
        % shape = cell(1,2);
        %%% CIRCLE
        % t = (0:2*pi/200:2*pi);
        % shape{1} = sin(t);
        % shape{2} = cos(t);
        %%% TRIANGLE (of same surface as circle)
        % lambda = 2*pi/(3*sqrt(3));
        % shape{1} = lambda*[-sqrt(3)/2,0,sqrt(3)/2];
        % shape{2} = lambda*[-0.5,1,-0.5];
        %%% SQUARE
        % shape{1} = [ -1 -1 1 1 ];
        % shape{2} = [ -1 1 1 -1 ];
        %
        % The figure can be saved with the following code:
        %
        % outfile = 'example.bmp';
        %%%%%% take a snapshot
        % f = getframe(gcf);
        % [img, map] = frame2im(f);
        % imwrite(img, strcat(outdir,outfile));
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% success is set to false if script encounters catastrophic
        %%% failure
        success = true;
        
        %%% specific maximum number of item locations
        maxnum = numerosity / totaloccupiedarea; %%% prepare a large enough number of locations
        
        %%%% verify item size information, and set it if necessary
        if (rmax==0)
            rmax=min( (1/sqrt(maxnum))/2.5 , 0.3); %%% this is the maximum radius of a given dot
            s=sprintf('Maximum item radius set to automatic value: %5.3f',rmax);
            disp(s);
        end
        
        %%%% precompute coordinates where dots can potentially appear
        % the coordinates are on a square matrix, but we select only those that fall
        % well within the circle
        % until we have enough locations to reach the maximum specified numerosity
        
        ndesired = floor(maxnum*1.3)+1; %%% larger number is used so that not all locations are filled by objects
        %%%% n is the side of a square matrix of possible locations;
        n=round(sqrt(ndesired));
        while 1
            coords = zeros((n+1)*(n+1),2);
            a=0;
            for i = 0:n
                for j=0:n
                    a=a+1;
                    coords(a,1)=(i-n/2 + (mod(j,2)-0.5)/4 ) / (n/2) ;
                    coords(a,2)=(j-n/2 + (mod(i,2)-0.5)/4 ) / (n/2);
                    %%% this generates a roughly square matrix, with offsets that
                    %%% prevent subjects from noticing the regular matrix
                end
            end
            dist = sqrt( coords(:,1).*coords(:,1) + coords(:,2).*coords(:,2) );
            coords = coords(dist<(0.98-1.2*rmax),:); %%% leave an empty ring of 0.02, and allow for some variability (20%) in central coordinates
            if (size(coords,1)>ndesired)
                break
            else
                n=n+1;
            end
        end
        % disp(n); %%% eventual number of locations chosen
        
        %%% sort the coordinates as a function of distance from center
        maxcoords = size(coords,1);
        dist = sqrt( coords(:,1).*coords(:,1) + coords(:,2).*coords(:,2) );
        [d,i]=sort(dist);
        coords = coords(i,:);
        
        %disp(sprintf('Side of square location matrix: %d',n));
        %disp(sprintf('Total number of locations generated: %d',maxcoords));
        
        if (0)
            %%% use this code to visualize the positions where items can appear
            hold off;
            fill(circle_x,circle_y,backcolor)	% large background circle
            %fill(square_x,square_y,backcolor)	% large background square
            axis square off image fill;
            hold on;
            r = 0.01;
            for i=1:maxcoords
                fill(circle_x*r+coords(i,1),circle_y*r+coords(i,2),forecolor);
            end
            pause;
        end
        %%% specify the ideal size and spacing parameters;
        
        %%% manipulation of the radius of dots
        
        r=sqrt(itemsize)*rmax;
        
        %%%% generate the stimulus
        
        usedcoords=max(round(totaloccupiedarea*maxcoords),1);
        if (usedcoords<numerosity) | (usedcoords>maxcoords)
            s=sprintf('Sorry, I cannot place %d objects amongst %d positions while occupying',numerosity,maxnum);
            disp(s);
            s=sprintf('only %5.1f percent of the display!',100*totaloccupiedarea);
            disp(s);
            disp('You can either increase the total area, or increase maxnum');
            error('NO DIPLAY GENERATED');
        end
        
        %%% scale up the display if necessary
        alpha = sqrt((max(usedcoords,3))/maxcoords);  %%% alpha<1, indicates ratio of occupied positions
        roffset= 0.9*rand*(1-alpha);
        % scalefact = 0.9*alpha;
        scalefact = 1;
        
        omega = 2*pi*rand;
        offset = [ roffset * sin(omega)  roffset * cos(omega) ];
        
        %%% shuffle to select a random subset of the used coords
        order = randperm(usedcoords);
        
        %%% add random rotation
        omega = 2*pi*rand;
        %        rotmat = [1 0 ; 0 1 ];
        rotmat = scalefact * [ sin(omega) cos(omega) ; - cos(omega) sin(omega) ];
        
        if useFill
            colordef none;      %set the color defaults to their MATLAB values
            backColor = [.5 .5 .5];
            ringColor = [1 1 1];
            centerColor = [0 0 0];
            
            pos=[200, 200, windowSize, windowSize]; %position of window, add 10 to size to crop later
            figure(101);
            clf;
            set(gcf,'Position',pos);
            set(gcf,'Color',backColor);
            axes('position',[0 0 1 1 ]);
            axis([-1 1 -1 1])
            % BACKGROUND CIRCLE
            t = (0:2*pi/200:2*pi);
            circle_shape{1} = sin(t);
            circle_shape{2} = cos(t);
            h = fill(circle_shape{1},circle_shape{2},backColor);	% large background circle
            set(h,'EdgeColor','none');
            %fill(square_x,square_y,backcolor)	% large background square
            axis square off image fill;
            
            % generate vectors of coordinates to draw the shapes
            if strcmp(shape,'circle') || strcmp(shape,'gabor');
                % CIRCLE
                t = (0:2*pi/200:2*pi);
                shapeCoords(:,1) = sin(t);
                shapeCoords(:,2) = cos(t);
            elseif strcmp(shape,'triangle');
                % TRIANGLE
                lambda = 2*pi/(3*sqrt(3)); % this is to make the surface of the triangle
                % equal to the surface of the circle
                shapeCoords(:,1) = lambda*[-sqrt(3)/2,0,sqrt(3)/2];
                shapeCoords(:,2) = lambda*[-0.5,1,-0.5];
            elseif strcmp(shape,'square');
                % SQUARE
                lambda = pi/4;
                shapeCoords(:,1) = lambda*[ -1 -1 1 1 ];
                shapeCoords(:,2) = lambda*[ -1 1 1 -1 ];
            else
                error('unknown shape %s',shape)
            end
            
            hold on
            for i=1:numerosity
                c = coords(order(i),:) * rotmat + ((rand(1,2)-0.5)*0.4)*r + offset;
                h(1) = fill(shapeCoords(:,1)*r+c(1),shapeCoords(:,2)*r+c(2),ringColor);
                set(h(1),'EdgeColor','none');
                r2 = sqrt(r^2/2);
                h(2) = fill(shapeCoords(:,1)*r2+c(1),shapeCoords(:,2)*r2+c(2),centerColor);
                set(h(2),'EdgeColor','none');
            end
            hold off
            f = getframe(gcf);
            img = frame2im(f);
            img = double(rgb2gray(img));
            img = (img)./range(img(:)); %scale to unit range
            img = img - mean(img(:)); %bring mean luminance to zero
            img = img/max(abs(img(:))); %Scale so max signed value is 1
            [~,minIdx] = min(abs([255,0]-lumIdx));
            if minIdx == 1
                img = (255-lumIdx)*img+lumIdx; % Scale into x-255 range
            else
                img = lumIdx*img+lumIdx; % Scale into 0-x range
            end
            img = uint8(img);
        else
            img = uint8(ones(windowSize, windowSize)*lumIdx);
            meshSize = windowSize;
            meshMapping = linspace(-1,1,meshSize);
            scaleFactor = 1;
            % radius of outer edge
            meshRadius(1) = scaleFactor*r.*meshSize/2;
            % radius of inner edge, corresponding to half the area
            meshRadius(2) = sqrt(meshRadius(1)^2/2);
            % round to whole numbers
            meshRadius = 2.*round(meshRadius./2);
            trueRadius = round(meshRadius./scaleFactor);
            % get total area of the square around the dot
            totalPix = (meshRadius(1)*2+1)^2;
            % compute rounding error and throw an error if more than 25%
            roundError = abs((pi*meshRadius(1).^2 - (2*pi*meshRadius(2).^2)))/(pi*meshRadius(1).^2);
            if roundError > .25
                error('rounding error more than 25%');
            else
            end
            for i=1:numerosity
                % get coordinates, following Dehaene's code
                c = coords(order(i),:) * rotmat + ((rand(1,2)-0.5)*0.4)*r + offset;
                % find closest actual x and y pixel coords
                meshCoords(1) = findClosest(c(1),meshMapping);
                meshCoords(2) = findClosest(c(2),meshMapping);
                [x,y] = meshgrid(1:(meshRadius(1)*2+1),1:(meshRadius(1)*2+1));
                if strcmp(shape,'circle')
                    if i == 1 % just generate the dot once
                        % first outer edge
                        tmpImg = ~(sqrt(((x-meshRadius(1)-1).^2 / meshRadius(1)^2) + ((y-meshRadius(1)-1).^2 / meshRadius(1)^2)) > 1);
                        % then add inner
                        tmpImg = tmpImg + ~(sqrt(((x-meshRadius(1)-1).^2 / meshRadius(2)^2) + ((y-meshRadius(1)-1).^2 / meshRadius(2)^2)) > 1);
                        
                        % compute central luminance values
                        bgPix = length(find(tmpImg == 0));
                        ringPix = length(find(tmpImg == 1));
                        centerPix = length(find(tmpImg == 2));
                        maxLum = 256;
                        minLum = -1;
                        while maxLum > 255
                            minLum = minLum+1;
                            maxLum = round(((lumIdx * (totalPix-bgPix)) - (centerPix * minLum) )./ringPix);   % solve for ring luminance, given bg = lumIdx and center = minLum
                        end
                        tmpImg = imresize(tmpImg,[trueRadius(1)*2+1,trueRadius(1)*2+1],'method','nearest');
                        tmpImg(tmpImg == 0) = lumIdx;
                        tmpImg(tmpImg == 1) = maxLum;
                        tmpImg(tmpImg == 2) = minLum;
                    else
                    end
                    curCoords = [(meshCoords(1)-trueRadius(1)):(meshCoords(1)+trueRadius(1)); (meshCoords(2)-trueRadius(1)):(meshCoords(2)+trueRadius(1))];
                    if any(curCoords(:) > windowSize) || any(curCoords(:) < 1) % take care of "spill-out"
                        success = false;
                        img = [];
                        return; 
                    else
                        imgChunk = img(curCoords(1,:),curCoords(2,:));
                        if ( unique(imgChunk(tmpImg == minLum)) == lumIdx ) || ( unique(imgChunk(tmpImg == maxLum)) == lumIdx ) % check for overlap
                            imgChunk(tmpImg == minLum) = minLum;
                            imgChunk(tmpImg == maxLum) = maxLum;
                            img(curCoords(1,:),curCoords(2,:)) = imgChunk;
                        else
                            success = false;
                            img = [];
                            return; 
                        end
                    end
                elseif strcmp(shape,'gabor');
                    gaborSize = (trueRadius(1)*2)+1;
                    gaborSigma = gaborSize/7;
                    tmpImg = gabor2d(gaborSize,gaborSize,gaborSize/360,rand*360,gaborSigma,gaborSigma);
                    tmpImg = (tmpImg)./range(tmpImg(:)); %scale to unit range
                    tmpImg = tmpImg - mean(tmpImg(:)); %bring mean luminance to zero
                    tmpImg = tmpImg/max(abs(tmpImg(:))); %Scale so max signed value is 1
                    [~,minIdx] = min(abs([255,0]-lumIdx));
                    if minIdx == 1
                        tmpImg = (255-lumIdx)*tmpImg+lumIdx; % Scale into x-255 range
                    else
                        tmpImg = lumIdx*tmpImg+lumIdx; % Scale into 0-x range
                    end
                    img((meshCoords(1)-trueRadius(1)):(meshCoords(1)+trueRadius(1)), (meshCoords(2)-trueRadius(1)):(meshCoords(2)+trueRadius(1))) = tmpImg;
                else
                end
                if abs(mean(tmpImg(:))-lumIdx) > 1 || max(tmpImg(:)) > 255 || min(tmpImg(:)) < 0
                    error('luminance matching failed during image generation');
                else
                end
            end
        end
    end

    function closestIdx = findClosest(val,list)
        [~,closestIdx] = min(abs(val-list));
    end

    function mat = gabor2d( sizeX, sizeY, freq, angle, sigmaX, sigmaY, meanX,meanY, dephase, cut)
        % gabor2d() - generate a two-dimensional gabor matrice.
        %
        % Usage:
        %   >> [ matrix ] = gabor2d(rows, columns);
        %   >> [ matrix ] = gabor2d( rows, columns, freq, ...
        %                             angle, sigmaR, sigmaC, meanR, meanC, dephase, cut)
        % Example :
        %	>> imagesc(gabor2d( 50, 50))
        %
        % Inputs:
        %   rows        - number of rows
        %   columns     - number of columns
        %   freq        - frequency of the sinusoidal function in degrees (default: 360/rows)
        %   angle       - angle of rotation of the resulting 2-D array in
        %                 degrees of angle {default: 0}.
        %   sigmaR      - standard deviation for rows {default: rows/5}
        %   sigmaC      - standard deviation for columns {default: columns/5}
        %   meanR       - mean for rows {default: center of the row}
        %   meanC       - mean for columns {default: center of the column}
        %   dephase     - phase offset in  degrees {default: 0}. A complex Gabor wavelet
        %                 can be build using gabor2dd(...., 0) + i*gabor2d(...., 90),
        %                 0 and 90 being the phase offset of the real and imaginary parts
        %   cut	        - percentage (0->1) of maximum value below which to remove values
        %                 from the matrix {default: 0}
        % Ouput:
        %   matrix - output gabor matrix
        %
        % Author: Arnaud Delorme, CNL / Salk Institute, 2001
        
        % Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
        %
        % This program is free software; you can redistribute it and/or modify
        % it under the terms of the GNU General Public License as published by
        % the Free Software Foundation; either version 2 of the License, or
        % (at your option) any later version.
        %
        % This program is distributed in the hope that it will be useful,
        % but WITHOUT ANY WARRANTY; without even the implied warranty of
        % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        % GNU General Public License for more details.
        %
        % You should have received a copy of the GNU General Public License
        % along with this program; if not, write to the Free Software
        % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
        
        if nargin < 2
            help gabor2d
            return;
        end;
        if nargin < 3
            freq = 360/sizeX;
        end;
        if nargin < 4
            angle = 0;
        end;
        if nargin < 5
            sigmaX = sizeX/5;
        end;
        if nargin < 6
            sigmaY = sizeY/5;
        end;
        if nargin < 7
            meanX = (sizeX+1)/2;
        end;
        if nargin < 8
            meanY = (sizeY+1)/2;
        end;
        if nargin < 9
            dephase = 0;
        end;
        if nargin < 10
            cut = 0;
        end;
        freq = freq/180*pi;
        
        X = linspace(1, sizeX, sizeX)'* ones(1,sizeY);
        Y = ones(1,sizeX)'   		  * linspace(1, sizeY, sizeY);
        %[-sizeX/2:sizeX/2]'*ones(1,sizeX+1);
        %Y = ones(1,sizeY+1)'   *[-sizeY/2:sizeY/2];
        
        rotatedmat = ((X-meanX)+i*(Y-meanY)) * exp(i*angle/180*pi);
        mat = sin(real(rotatedmat)*freq + dephase/180*pi).*exp(-0.5*(  ((X-meanX)/sigmaX).*((X-meanX)/sigmaX)...
            +((Y-meanY)/sigmaY).*((Y-meanY)/sigmaY)))...
            /((sigmaX*sigmaY)^(0.5)*pi);
        
        if cut > 0
            maximun = max(max(mat))*cut;
            I = find(mat < maximun);
            mat(I) = 0;
        end;
    end

    function TestSubFunction
        disp( 'inside TestSubfunction' )
    end

    function rV = GetParamArray( aPartName, aParamName )
        % *** Sin Step Types not included in logic below? ***
        
        % For the given part and parameter name, return an array of values
        % corresponding to the steps in a sweep.  If the requested param is
        % not swept, the array will contain all the same values.
        
        % tSpatFreqSweepValues = GetParamArray( '1', 'Spat Freq (cpd)' );
        
        % Here's an example of sweep type specs...
        %
        % definitions{end-2} =
        % 	{
        % 		'Fixed'         'constant'   { }
        % 		'Contrast'      'increasing' { { '1' 'Contrast (pct)' } { '2' 'Contrast (pct)' } }
        % 		'Spat Freq'      'increasing' { { '1' 'Spat Freq (cpd)' } { '2' 'Spat Freq (cpd)' } }
        % 	}
        
        tNCStps    = xDiva.GrabCellValue( timing, 'nmbCoreSteps' );
        tSweepType = xDiva.GrabCellValue( parameters{iS}, 'Sweep Type' );
        
        % we need to construct a swept array if any of the {name,value} in definitions{iSweep}{:,3}
        
        sweepList = xDiva.GrabCellValue( definitions{iSweep}, tSweepType, 3 );		% { {part,param}, {part,param} ... }
        
        % check for sweep
        % determine if any definitions{iSweep}{ iRow, { {part,param}... } } match arguments aPartName, aParamName
        partMatch  = ismember( aPartName,  cellfun( @(x)x{1}, sweepList, 'UniformOutput', false ) ); % will be false for "'Fixed' 'constant' {}"
        paramMatch = ismember( aParamName, cellfun( @(x)x{2}, sweepList, 'UniformOutput', false ) );
        if partMatch && paramMatch
            tSweepStart = xDiva.GrabCellValue( parameters{iS}, 'Sweep Start' );
            tSweepEnd   = xDiva.GrabCellValue( parameters{iS}, 'Sweep End' );
            if strcmpi( xDiva.GrabCellValue( parameters{iS}, 'Step Type' ), 'Lin Stair' );
                rV = linspace( tSweepStart, tSweepEnd, tNCStps )';
            else
                rV = logspace( log10(tSweepStart), log10(tSweepEnd), tNCStps )';
            end
        else
            rV = repmat( xDiva.GrabCellValue( parameters{eval(['i',aPartName])}, aParamName ), [ tNCStps, 1 ] );
        end
        
    end


end











