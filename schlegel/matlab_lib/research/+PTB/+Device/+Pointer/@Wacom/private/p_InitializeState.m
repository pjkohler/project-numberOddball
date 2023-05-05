function p_InitializeState(wac)
% p_InitializeState
% 
% Description:	initialize the state of the wacom devices
% 
% Syntax:	p_InitializeState(wac)
% 
% Updated: 2012-07-20
% Copyright 2012 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
[h,sz,rect,szVA] = wac.parent.Window.Get('main');

[xStylus,yStylus,butStylus,pStylus,txStylus,tyStylus]	= p_GetStylus(wac,sz);
[xEraser,yEraser,butEraser,pEraser,txEraser,tyEraser]	= p_GetEraser(wac,sz);
[xTouch,yTouch]											= p_GetTouch(wac,sz);

wac.last.stylus	= struct(...
					'x'			, xStylus	, ...
					'y'			, yStylus	, ...
					'button'	, butStylus	, ...
					'p'			, pStylus	, ...
					'tx'		, txStylus	, ...
					'ty'		, tyStylus	  ...
					);
wac.last.eraser	= struct(...
					'x'			, xEraser	, ...
					'y'			, yEraser	, ...
					'button'	, butEraser	, ...
					'p'			, pEraser	, ...
					'tx'		, txEraser	, ...
					'ty'		, tyEraser	  ...
					);
wac.last.touch	= struct(...
					'x'			, xTouch	, ...
					'y'			, yTouch	, ...
					'button'	, 0			, ...
					'p'			, 0			, ...
					'tx'		, 0			, ...
					'ty'		, 0			  ...
					);
