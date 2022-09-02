function [psdiff]=ps_diff (inp,theta)

[~,~,PS1,~] = lklfn2(inp,theta);
JN=inp.JN;
psdiff=JN-PS1;
end
