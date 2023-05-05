function [b,cPathTract] = FSLProbtrackx(cDirSubject,cPathSeed,varargin)
% FSLProbtrackx
% 
% Description:	call probtrackx. see FSL's probtrackx help for details about the
%				optional parameters.
% 
% Syntax:	[bSuccess,cPathTract] = FSLProbtrackx(cDirSubject,cPathSeed,<options>)
% 
% In:
%	cDirSubject	- the subject's FSL DTI directory, or a cell of directories to
%				  run probtracx for multiple subjects.
%	cPathSeed	- the path to a binary NIfTI seed volume for single mask
%				  tractography, a cell of paths for multiple mask tractography,
%				  or a cell of cells of paths for multiple probtrackx runs.
%				  seeds must be in diffusion space. can also be label files if
%				  you are following the procedure here:
%				  http://fsl.fmrib.ox.ac.uk/fsl/fsl-4.1.9/fdt/fdt_surface.html.
%	<options>:
%		name:				(<auto>) the name to use for the output directory,
%							or a cell of names
%		dir_out:			(<input>.probtrackX/<name>/) the output directory
%							path, or a cell of paths. overrides <name>.
%		mesh:				(<none>) for label-based tractography, the path/cell
%							of paths to the ascii-formatted surface mesh on
%							which the labels are based
%		seedref:			(<none>) the path/cell of paths to the reference
%							volume for the seeds. required if seeds are not in
%							diffusion space.
%		xfm:				(<none>) the path/cell of paths to FLIRT/FNIRT
%							transforms from seed space to diffusion space.
%							required if seeds are not in diffusion space.
%		invxfm:				(<none>) the path/cell of paths to FLIRT/FNIRT
%							inverse transforms from diffusion space to seed
%							space. required if seeds are not in diffusion space
%							and FNIRT transforms are used.
%		waypoint:			(<none>) the path/cell of paths/cell of cells of
%							paths to binary NIfTI volumes to use as waypoint
%							masks (i.e. tracts must pass through all of these
%							masks to be kept)
%		exclusion:			(<none>) the path/cell of paths/cell of cells of
%							paths to binary NIfTI volumes to use as exclusion
%							masks (i.e. any tract passing through any of these
%							masks will be rejected). note that these are ORed
%							into a single mask saved in the output directory
%							since probtrackx only takes one exclusion mask.
%		termination:		(<none>) the path/cell of paths/cell of cells of
%							paths to binary NIfTI volumes to use as termination
%							masks (i.e. any tract hitting any of these masks
%							will be stopped). note that these are ORed into a
%							single mask saved in the output directory since
%							probtrackx only takes one termination mask.
%		classification:		(<none>) the path/cell of paths/cell of cells of
%							paths to binary NIfTI volumes to use as
%							classification masks (i.e. seeds_to_<name> volumes
%							are created for each mask specified here). this only
%							works if a single seed is specified.
%		lengthcorrect:		(false) true to multiply each voxel's output value
%							by the expected length of tracts that cross it.
%							somewhere on the FSL mailing list it is mentioned
%							that this isn't recommended for quantitative
%							studies, only for classification. (note this is
%							probtrackx's pd option).
%		nsample:			(5000) the number of samples to send out from each
%							voxel in each seed
%		nstep:				(2000) the number of steps per sample
%		threshdistance:		(0) reject tracts that are below this length, in mm
%		threshcurvature:	(0.2) the curvature threshold for choosing each
%							sample's next position at each step
%		threshfiber:		(0.01) probtrackx's fibthresh option (not sure what
%							this is)
%		steplength:			(0.5) the length of each sample step, in mm
%		loopcheck:			(true) true to check for and terminate looping paths
%		usef:				(false) use anisotropy to constrain tracking
%		modeuler:			(false) true to use modified euler streaming. this
%							is apparently more accurate but takes longer.
%		rseed:				(<none>) the random seed value for tracking
%		cores:				(1) the number of processor cores to use
%		force:				(true) true to force probtrackx to run even if the
%							output files already exist
%		force_pre:			(false) true to force preprocessing (e.g.
%							calculation of transforms, etc.)
%		silent:				(false) true to suppress status messages
% 
% Out:
% 	bSuccess		- true if the process ran successfully
%	strNameTract	- the tract name
%	strScript		- the script for the call to probtrackx
% 
% Updated: 2015-05-01
% Copyright 2015 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
bSuccess	= false;

opt	= ParseArgs(varargin,...
		'name'				, []			, ...
		'dirout'			, []			, ...
		'seedspace'			, 'diffusion'	, ...
		'ref'				, []			, ...
		'xfm'				, []			, ...
		'xfminv'			, []			, ...
		'wmstopmask'		, true			, ...
		'wm_grow'			, -1			, ...
		'dir_fs'			, []			, ...
		'waypoint'			, {}			, ...
		'exclusion'			, {}			, ...
		'termination'		, {}			, ...
		'classification'	, {}			, ...
		'nsample'			, 5000			, ...
		'nstep'				, 2000			, ...
		'steplength'		, 0.5			, ...
		'threshcurvature'	, 0.2			, ...
		'threshdistance'	, 0				, ...
		'lengthcorrect'		, false			, ...
		'loopcheck'			, true			, ...
		'usef'				, false			, ...
		'modeuler'			, false			, ...
		'rseed'				, []			, ...
		'force'				, true			, ...
		'forceprep'			, false			, ...
		'dotract'			, true			, ...
		'silent'			, false			  ...
		);

strDirSubject	= AddSlash(strDirSubject);

cPathSeed	= ForceCell(cPathSeed);
nSeed		= numel(cPathSeed);

%output files
	if isempty(opt.dirout)
		strDirPTX	= FSLDirPTX(strDirSubject);
		
		if isempty(opt.name)
			[strDirOut,strNameTract]	= GetUniqueDir(strDirPTX,'create',false);
		else
			strNameTract	= opt.name;
			strDirOut		= DirAppend(strDirPTX,strNameTract);
		end
	else
		strNameTract	= char(DirSplit(opt.dirout,'limit',1));
		strDirOut		= opt.dirout;
	end
	
	CreateDirPath(strDirOut,'error',true);
	
	strPathTract			= PathUnsplit(strDirOut,'fdt_paths','nii.gz');
	strPathWaytotalFirst	= deal(PathUnsplit(strDirOut,'waytotal'));
	
	%add suffixes depending on options
		strSuffix						= FSLSuffixTract('seedspace',opt.seedspace,'lengthcorrect',opt.lengthcorrect);
		[strPathTract,strPathWaytotal]	= varfun(@(f) PathAddSuffix(f,strSuffix,'favor','nii.gz'),strPathTract,strPathWaytotalFirst);
	
	if opt.dotract && ~opt.force && all(FileExists({strPathTract strPathWaytotal}))
		bSuccess	= true;
		return;
	end

strDirBedpostx	= [RemoveSlash(strDirSubject) '.bedpostX'];
strBaseSample	= PathUnsplit(strDirBedpostx,'merged');
strPathNoDif	= PathUnsplit(strDirSubject,'nodif_brain_mask','nii.gz');

[opt.waypoint,opt.exclusion,opt.termination,opt.classification]	= ForceCell(opt.waypoint,opt.exclusion,opt.termination,opt.classification);

%get the transform info
	bTransform	= true;
	
	[strPathRef,strPathXFM,strPathXFMInv]	= deal([]);
	
	switch lower(opt.seedspace)
		case 'freesurfer'
		%calculate the freesurfer->fa transform
			strDirMRI	= DirAppend(opt.dir_fs,'mri');
			strPathRef	= PathUnsplit(strDirMRI,'brain','nii.gz');
			
			[b,strPathXFM,strPathXFMInv]	= FreeSurfer2FA(opt.dir_fs,strDirSubject,'force',opt.forceprep,'silent',opt.silent);
			if ~b
				status('Could not calculate the FreeSurfer<-->FA registration.','warning',true,'silent',opt.silent); 
				return;
			end
		case 'diffusion'
			bTransform	= false;
		otherwise
	end
	
	%keep the user specified values
		strPathRef		= unless(opt.ref,strPathRef);
		strPathXFM		= unless(opt.xfm,strPathXFM);
		strPathXFMInv	= unless(opt.xfminv,strPathXFMInv);
%save the inverse white matter mask
	if opt.wmstopmask
		[b,strPathWMInvFS]	= FreeSurferMaskWMInverse(opt.dir_fs,...
								'grow'		, opt.wm_grow	, ...
								'force'		, opt.forceprep	, ...
								'silent'	, opt.silent	  ...
								);
		if ~b
			return;
		end
		
		%transform the mask and copy it to the subject's FSL directory
			strPathWMInv	= PathAddSuffix(PathUnsplit(strDirSubject,PathGetFileName(strPathWMInvFS)),['-' lower(opt.seedspace)],'favor','nii.gz');
			
			if opt.forceprep || ~FileExists(strPathWMInv)
				if isempty(opt.seedspace)
					error('The seed space must be defined if the white matter stop mask is used.');
				end
				
				switch lower(opt.seedspace)
					case 'freesurfer'
					%copy
						if ~copyfile(strPathWMInvFS,strPathWMInv) && ~FileExists(strPathWMInv)
							status('Could not copy the inverse white matter mask to the output directory','warning',true,'silent',opt.silent);
							return;
						end
					case 'diffusion'
						%calculate the freesurfer->fa transform
							[b,strPathXFM]	= FreeSurfer2FA(opt.dir_fs,strDirSubject,'force',opt.forceprep,'silent',opt.silent);
							if ~b
								status('Could not calculate the FreeSurfer<-->FA registration.','warning',true,'silent',opt.silent); 
								return;
							end
						%transform the white matter mask
							strPathRef	= PathUnsplit(strDirSubject,'nodif_brain_mask','nii.gz');
							
							if ~FSLXFM(strPathWMInvFS,strPathXFM,strPathRef,'output',strPathWMInv,'mask',true,'force',opt.forceprep,'silent',opt.silent)
								status('Could not transform the inverse white matter mask to FA space.','warning',true,'silent',opt.silent); 
								return;
							end
					otherwise
						error('The seed space must be either ''freesurfer'' or ''diffusion'' if the white matter stop mask is used.');
				end
			end
		
		opt.termination{end+1}	= strPathWMInv;
	end
%save the seed, waypoint, and classification mask lists
	bSeed		= nSeed>0;
	bNetwork	= nSeed>1;
	if bNetwork
		strPathSeed	= PathUnsplit(strDirOut,'seed','txt');
		%probtrackx doesn't work unless there's a trailing line feed.  guess how
		%long it took me to figure this out...
		fput([join(cPathSeed,10) 10],strPathSeed);
	elseif bSeed
		strPathSeed	= cPathSeed{1};
	else
		strPathSeed	= '';
	end
	
	strPathWaypoint	= PathUnsplit(strDirOut,'waypoint','txt');
	bWaypoint		= ~isempty(opt.waypoint);
	if bWaypoint
		fput([join(opt.waypoint,10) 10],strPathWaypoint);
	end
	
	strPathClassification	= PathUnsplit(strDirOut,'classification','txt');
	bClassification			= ~isempty(opt.classification);
	if bClassification
		if nSeed>1
			error('Classification masks can only be used if one seed mask is specified.');
		end
		
		fput([join(opt.classification,10) 10],strPathClassification);
	end
%save the ORed other exclusion/termination masks
	%exclusion
		strPathExclusion	= PathUnsplit(strDirOut,'exclusion','nii.gz');
		bExclusion			= ~isempty(opt.exclusion);
		if bExclusion && ~MRIMaskMerge(opt.exclusion,strPathExclusion,'force',opt.force,'silent',opt.silent)
			status('Could not OR exclusion masks','warning',true,'silent',opt.silent);
			return;
		end
	%termination
		strPathTermination	= PathUnsplit(strDirOut,'termination','nii.gz');
		bTermination		= ~isempty(opt.termination);
		if bTermination && ~MRIMaskMerge(opt.termination,strPathTermination,'force',opt.force,'silent',opt.silent)
			status('Could not OR termination masks','warning',true,'silent',opt.silent);
			return;
		end
%run probtrackx
	bRSeed			= ~isempty(opt.rseed);

	strOut			= PathGetFilePre(strPathTract,'favor','nii.gz');
	strDistThresh	= num2str(opt.threshdistance);
	strStepLength	= num2str(opt.steplength);
	strRSeed		= num2str(opt.rseed);
	
	if bSeed
		%temporarily rename the existing waytotal file
			bWTExist	= opt.dotract && FileExists(strPathWaytotalFirst);
			if bWTExist
				strPathWaytotalTemp	= PathAddSuffix(strPathWaytotalFirst,'_previous');
				movefile(strPathWaytotalFirst,strPathWaytotalTemp);
			end
		
		b	= CallProcess('probtrackx',[{}...
				'-x'	, strPathSeed				, ...
				'-s'	, strBaseSample				, ...
				'-o'	, strOut					, ...
				'-m'	, strPathNoDif				, ...
				'-P'	, opt.nsample				, ...
				'-S'	, opt.nstep					, ...
				'-c'	, opt.threshcurvature		, ...
				['--dir="' strDirOut '"']			, ...
				['--distthresh=' strDistThresh]	, ...
				['--steplength=' strStepLength]	, ...
				'--forcedir'						, ...
				'--opd'								, ...
				'--mode=seedmask'					, ...
				conditional(bNetwork			, '--network'										, []) , ...
				conditional(opt.lengthcorrect	, '--pd'											, []) , ...
				conditional(bClassification	, '--os2t'											, []) , ...
				conditional(opt.loopcheck		, '-l'												, []) , ...
				conditional(opt.usef			, '-f'												, []) , ...
				conditional(opt.modeuler		, '--modeuler'										, []) , ...
				conditional(bTransform			, ['--seedref="' strPathRef '"']					, []) , ...
				conditional(bTransform			, ['--xfm="' strPathXFM '"']						, []) , ...
				conditional(bTransform			, ['--invxfm="' strPathXFMInv '"']					, []) , ...
				conditional(bWaypoint			, ['--waypoints="' strPathWaypoint '"']			, []) , ...
				conditional(bTermination		, ['--stop="' strPathTermination '"']				, []) , ...
				conditional(bClassification	, ['--targetmasks="' strPathClassification '"']	, []) , ...
				conditional(bExclusion			, ['--avoid="' strPathExclusion '"']				, []) , ...
				conditional(bRSeed				, ['--rseed=' strRSeed]							, []) ],...
				'run'		, opt.dotract	, ...
				'silent'	, opt.silent	  ...
				);
		
		if opt.dotract
			b	= ~b;
			
			if ~isequal(strPathWaytotal,strPathWaytotalFirst)
			%rename the waytotal file
				movefile(strPathWaytotalFirst,strPathWaytotal);
				
				if bWTExist
				%move the old waytotal back
					movefile(strPathWaytotalTemp,strPathWaytotalFirst);
				end
			end
		else
			bSuccess	= b{1};
			return;
		end
		
		if ~b
			status('Probtrackx failed','warning',true,'silent',opt.silent);
			return;
		end
	end
%success!
	bSuccess	= true;
