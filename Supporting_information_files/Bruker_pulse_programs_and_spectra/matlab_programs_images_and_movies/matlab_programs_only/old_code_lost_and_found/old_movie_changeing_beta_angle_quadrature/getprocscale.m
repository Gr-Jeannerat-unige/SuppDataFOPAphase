function [ nc]=getprocscale(path)
str_in=[path  '/procs'];
    %read acqus
    content=textread(str_in,'%s' );
    %sw         = str2num(char(content(strmatch('##$SW=',content')+1)));
    %pulprog    = char(content(strmatch('##$PULPROG=',content')+1));pulprog = pulprog(2:(end-1));
    % this is the number of point in F2, (2 td)
    NC_proc         = str2num(char(content(strmatch('##$NC_proc=',content')+1)));
  NC_diff=-14-NC_proc;nc=power(2,NC_diff);
  