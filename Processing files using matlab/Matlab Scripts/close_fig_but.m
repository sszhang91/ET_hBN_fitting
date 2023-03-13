function close_fig_but(figs2keep)
    %takes as input figs2keep as vector of numbers
    %eg [1, 3, 4]
    
    % Uncomment the following to 
    % include ALL windows, including those with hidden handles (e.g. GUIs)
    % all_figs = findall(0, 'type', 'figure');
    all_figs=findobj(0,'type','figure');
    if nargin>0
        delete(setdiff(all_figs,figs2keep))
    else
        delete(all_figs)
    end
end