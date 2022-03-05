function out = TDR(Xin, Yin, conds, varargin)
% Input:
%   Xin: regressors encoding task variables(t.v.). 
%        Size of #trial X #t.v.
%        (e.g., 1 task variable (e.g., amount of reward) can be encoded as
%         series of ones and zeros. 1 for high reward and 0 for low reward.)
%        NOTES. scale of each task variable should be (ideally) same to avoid some
%        variable specific "leaning" of regression coefficients.
%   Yin: neural firing rates. 
%        Size of #time X  #trial or #neuron X #time X #trial
%   varargin:
%       'normY' or 'normy' : zscoring of each neuron. (defaults: off)
%                            (mean and std are computed after combining
%                             time and trials.)
%       'add_Intcpt'       : no use of intercept in the GLM process. (defaults: on)
%       'add_Interact'     : add interaction term.

do_zscr = false;
add_Intcpt = true; 
add_Interact = off;

for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'normY' | 'normy'
                % zscoring neural firing rates.
                do_zscr = true;
            case 'noIntcpt'
                add_Intcpt = false;
        end
    end
end


if numel(size(Yin)) == 3
    
elseif numel(size(Yin)) == 2
    
else
    error('Inadequate size of neural firing rates!');
end







end