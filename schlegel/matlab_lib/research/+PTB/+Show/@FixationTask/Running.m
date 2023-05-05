function b = Running(ft)
% PTB.FixationTask.Running
% 
% Description:	test whether the fixation task is running
% 
% Syntax:	b = ft.Running
%
% Out:
%	b	- true if the fixation task is running
% 
% Updated: 2011-12-20
% Copyright 2011 Alex Schlegel (schlegel@gmail.com).  This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported
% License.
b	= ft.parent.Scheduler.Running('fixation_task');
